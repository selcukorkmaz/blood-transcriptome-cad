# =========================================
# Validation Phase: GSE113079 Differential Expression
# =========================================

# Set working directory
setwd("~/GitHub/cad")

# ─── 1. Load libraries ───────────────────────────────────────────
library(GEOquery)
library(data.table)
library(dplyr)           # For data wrangling
library(stringr)         # For handling string patterns if needed
library(biomaRt)
library(ggplot2)
library(limma)
library(patchwork)
library(EnhancedVolcano)
library(umap)

# ─── 2. Download and load GEO dataset ───────────────────────────
gse113079_list <- getGEO("GSE113079", GSEMatrix = TRUE)

# Typically, you get one main ExpressionSet for a single-platform study:
length(gse113079_list)
# [1] 1 (if there's only one, we'll use gse113079_list[[1]])

gse113079 <- gse113079_list[[1]]

# View sample (patient) metadata
# pData(gse113079) %>% head()

# ─── 3. Annotate sample metadata and define groups ───────────────
meta <- pData(gse113079)
meta$group <- ifelse(
  grepl("CAD patient", meta$title),  "CAD",
  ifelse(grepl("Healthy control", meta$title), "Control", NA)
)
# Check distribution:
table(meta$group)
pData(gse113079) <- meta

# ─── 4. Extract expression matrix and QC ────────────────────────
expr_matrix <- exprs(gse113079)
dim(expr_matrix)
# e.g. 100k probes x 141 samples

# Ensure columns match metadata rows
all(colnames(expr_matrix) == rownames(meta))  # Should be TRUE

# Quick check & skip any log2
qx   <- as.numeric(quantile(expr_matrix, c(0,0.25,0.5,0.75,0.99,1)))
logC <- (qx[5] > 100) ||
  ((qx[6] - qx[1] > 50) && (qx[2] > 0))
if (logC) {
  stop("These look like raw intensities — you need to log2(x+1) first!")
} else {
  message("Data already on log2 scale, skipping log transform.")
}

gse113079_exprs <- as.data.frame(expr_matrix %>% t())
gse113079_exprs$group <- meta$group

write.table(gse113079_exprs, file="data/gse113079_val_expr.txt",    sep="\t", quote=FALSE, row.names=FALSE)

meta$Group <- factor(meta$group, levels = c("Control", "CAD"))

pca_obj <- prcomp(t(expr_matrix), scale. = TRUE)
pca_df <- as.data.frame(pca_obj$x[, 1:2])          # PC1 & PC2 scores
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Group <- meta$Group
summary(pca_obj)


# coloured by biology
ggplot(pca_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 3) +
  labs(title = "Validation PCA – coloured by group") +
  theme_minimal()
ggsave("Figures/Val_PCA_group.png", width = 6, height = 5)


set.seed(123)                                  # reproducible layout
ump <- umap(t(expr_matrix))                    # use same matrix as limma

umap_df <- data.frame(
  UMAP1 = ump$layout[, 1],
  UMAP2 = ump$layout[, 2],
  Group = meta$Group
)
umap_df

ggplot(umap_df, aes(UMAP1, UMAP2, colour = Group)) +
  geom_point(size = 3) +
  labs(title = "Validation UMAP – coloured by group") +
  theme_minimal()
ggsave("Figures/Val_UMAP_group.png", width = 6, height = 5)

# ─── 6. Build design matrix ────────────────────────────────────
group  <- factor(meta$group, levels = c("Control", "CAD"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)  # "Control", "CAD"

# ─── 7. Fit linear model and compute DE ─────────────────────────
fit     <- lmFit(expr_matrix, design)
cont.matrix <- makeContrasts(
  CADvsControl = CAD - Control,
  levels       = design
)
fit2    <- contrasts.fit(fit, cont.matrix)
fit2    <- eBayes(fit2)

# Get top differentially expressed probes
results <- topTable(
  fit2, coef = "CADvsControl",
  number = Inf, adjust.method = "fdr",
  sort.by = "p"
)
results %>% head()

# ─── 8. Annotate probes with featureData ───────────────────────
annot    <- fData(gse113079)
head(annot)

annot_df <- annot %>%
  dplyr::select(PROBE_ID = 1, GeneSymbol, everything()) %>%  # adapt as needed
  mutate(PROBE_ID = rownames(annot))

results$PROBE_ID    <- rownames(results)
results_annot       <- left_join(results, annot_df, by = "PROBE_ID")
head(results_annot)

# ─── 9. Load discovery-phase hits for merging ───────────────────
gse208194_hits <- read.table(
  "data/degs_annotated_gse208194.txt",
  header = TRUE, sep = "\t"
)
colnames(gse208194_hits)[8] <- "GeneSymbol"
head(gse208194_hits)

# ─── 10. Subset lncRNAs and mRNAs in discovery hits ─────────────
lnc_df  <- subset(
  gse208194_hits,
  grepl("lncRNA", gene_biotype) |
    gene_biotype %in% c("lincRNA", "antisense", "processed_transcript")
)
head(lnc_df)

mrna_df <- subset(
  gse208194_hits,
  gene_biotype == "protein_coding"
)
head(mrna_df)

# ─── 11. Prepare results_annot for merging ──────────────────────
results_annot_new <- results_annot[, c(
  "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
  "PROBE_ID", "GeneSymbol", "SPOT_ID", "GeneName",
  "EnsemblID", "lncRNA_ProbeName", "lncRNA ID", "gene"
)]
head(results_annot_new)
head(results_annot_new[, c("GeneSymbol", "gene")], 20)

results_annot_new$gene_ids <- ifelse(
  grepl("^ENSG", results_annot_new$gene),
  sub("\\..*", "", results_annot_new$gene),
  results_annot_new$gene
)
head(results_annot_new)

# ─── 12. Query Ensembl for HGNC symbols ─────────────────────────
library(biomaRt)
mart     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensg_ids <- results_annot_new$gene_ids

genes_ensbl <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensg_ids,
  mart       = mart
)
colnames(genes_ensbl) <- c("gene_ids", "GeneSymbol")
head(genes_ensbl)
genes_ensbl[genes_ensbl$GeneSymbol == "MIR663AHG", ]

# ─── 13. Query Ensembl by chromosomal region ───────────────────
results_annot$GenomicCoordinates <- sub("^chr", "", results_annot$GenomicCoordinates)

genes_coordinates <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "chromosomal_region",
  values     = results_annot$GenomicCoordinates,
  mart       = mart
)
colnames(genes_coordinates) <- c("gene_ids", "GeneSymbol")
head(genes_coordinates)

# ─── 14. Combine symbol sources and merge ──────────────────────
genes_combined <- rbind(genes_ensbl, genes_coordinates)
head(genes_combined)
head(results_annot_new)

merged_results <- left_join(
  results_annot_new, genes_combined,
  by     = "gene_ids",
  suffix = c("", "_new")
) %>%
  mutate(
    GeneSymbol = ifelse(
      is.na(GeneSymbol) | GeneSymbol == "",
      GeneSymbol_new,
      GeneSymbol
    )
  ) %>%
  dplyr::select(-GeneSymbol_new) %>%
  filter(!is.na(GeneSymbol))
# Final dataset in 'merged_results'

genes_combined[genes_combined$gene_ids == "ENSG00000227195", ]
results_annot_new[results_annot_new$gene_ids == "ENSG00000227195", ]
merged_results[merged_results$GeneSymbol == "MIR663AHG", ]

# ─── 15. Merge with discovery-phase results and find common hits ─
results_annot_new_merged <- merge(
  results_annot_new,
  genes_combined,
  by       = "gene_ids",
  all.x    = TRUE,
  suffixes = c("", "_updated")
)

common_hits <- inner_join(
  gse208194_hits, merged_results,
  by      = "GeneSymbol",
  suffix  = c(".discov", ".val")
)

best_probes <- common_hits %>%
  group_by(GeneSymbol) %>%
  slice_min(order_by = P.Value.val, with_ties = FALSE) %>%
  ungroup()


best_probes <- best_probes %>%
  filter(adj.P.Val.val < 0.05,
         logFC.discov * logFC.val > 0)

# head(best_probes)
# head(best_probes)
# dim(best_probes)
# table(best_probes$gene_biotype)
# colnames(best_probes)
# best_probes$GeneSymbol
# best_probes %>% print(n = 1000)
# best_probes[, c("GeneSymbol", "gene_biotype")] %>% print(n = 500)
# 
# # ─── 16. Pathway enrichment on validated mRNAs ──────────────────
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# geneList <- best_probes$GeneSymbol[best_probes$gene_biotype != "lncRNA"]
# genes_entrez <- bitr(
#   geneList,
#   fromType = "SYMBOL",
#   toType   = "ENTREZID",
#   OrgDb    = org.Hs.eg.db
# )
# # head(genes_entrez)
# 
# go_results <- enrichGO(
#   gene          = genes_entrez$ENTREZID,
#   OrgDb         = org.Hs.eg.db,
#   ont           = "ALL",
#   pAdjustMethod = "BH",
#   pvalueCutoff  = 0.05,
#   qvalueCutoff  = 0.05,
#   readable      = TRUE
# )
# # head(go_results)
# 
# kegg_results <- enrichKEGG(
#   gene         = genes_entrez$ENTREZID,
#   organism     = "hsa",
#   pAdjustMethod= "BH",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# )
# kegg_results <- setReadable(
#   kegg_results,
#   OrgDb   = org.Hs.eg.db,
#   keyType = "ENTREZID"
# )
# # head(kegg_results)
# 
# # Dotplots
# enrichplot::dotplot(go_results,   showCategory = 10) + ggtitle("GO Enrichment Analysis")
# enrichplot::dotplot(kegg_results, showCategory = 10) + ggtitle("KEGG Pathway Analysis")
# 
# library(enrichplot)
# emapplot(pairwise_termsim(go_results))
# emapplot(pairwise_termsim(kegg_results))
# 
# # Final check
# head(best_probes)

validated_genes = best_probes[,c("GeneSymbol", "gene_biotype", "logFC.discov", "adj.P.Val.discov", "logFC.val", "adj.P.Val.val")]
validated_genes %>% print(n=10000)

write.table(best_probes, "data/best_probes.txt", quote = F, row.names = F, sep = "\t")
write.table(validated_genes, "data/validated_genes.txt", quote = F, row.names = F, sep = "\t")


validation_set <- best_probes[,c("GeneSymbol", "gene_biotype", "logFC.discov", "adj.P.Val.discov", "logFC.val", "adj.P.Val.val")]
write.table(validation_set, "data/validation_set.txt", quote = F, row.names = F, sep = "\t")
