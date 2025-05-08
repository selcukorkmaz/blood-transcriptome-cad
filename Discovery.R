# =========================================
# Discovery Phase: Differential Expression
#    TPM data via limma (no quantile norm)
# =========================================

# Set working directory
setwd("~/GitHub/cad")

# ─── Load libraries ───────────────────────
library(data.table)      # Fast data I/O
library(dplyr)           # Data wrangling
library(biomaRt)         # Gene ID mapping
library(GEOquery)        # Fetch metadata
library(limma)           # Differential expression
library(clusterProfiler) # Enrichment analysis
library(org.Hs.eg.db)    # Human gene annotation
library(EnhancedVolcano) # Volcano plots
library(ggplot2)         # Ploresing
library(umap)            # UMAP visualization
library(ggfortify)     # for autoplot(); install if needed
library(patchwork)


# ─── 1. Read & subset TPM matrix ──────────
tpm_data <- fread("GSE208194_RawTPM.csv.gz", data.table = FALSE)
rownames(tpm_data) <- tpm_data[[1]]; tpm_data[[1]] <- NULL

keep_samps <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10",
                "S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                "S21","S22","S23","S24","S26","S27","S28","S29","S30","S31",
                "S32","S33","S34","S35","S36","S37","S38","S39","S40","S41",
                "S42","S44","S45","S46","S47","S48","S49","S50","S51","S52",
                "S53","S54","S55","S56","S57","S58","S59","S60","S95","S61",
                "S62","S63","S64","S65","S66","S67","S68","S69","S70","S71",
                "S72","S73","S74","S75","S77","S78","S79","S80","S81","S82",
                "S83","S84","S85","S86","S87","S88","S89","S90","S96")

tpm_data <- tpm_data[, keep_samps]

# ─── 2. Sample metadata & define groups ──────────
gse <- getGEO("GSE208194", GSEMatrix = TRUE)
meta <- pData(gse[[1]])

meta = meta %>% 
  filter(description %in% colnames(tpm_data))

groupInfo <- meta %>% 
  dplyr::select(description, raw = `treatment:ch1`) %>%
  mutate(group = case_when(
    raw %in% c("HF-","HF+") ~ "CAD",
    raw == "HVOL"          ~ "Control",
    TRUE                   ~ NA_character_
  )) 

# Confirm sample order matches exactly
tpm_data <- tpm_data[, groupInfo$description]

group <- factor(groupInfo$group, levels = c("Control", "CAD"))

# ─── 3. Map Ensembl IDs → gene symbols ───────────
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids_clean <- sub("\\..*", "", rownames(tpm_data))
gene_info <- getBM(
  attributes  = c("ensembl_gene_id", "external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids_clean,
  mart       = ensembl
)

rownames(tpm_data) <- ensembl_ids_clean
keep_ids <- intersect(rownames(tpm_data), gene_info$ensembl_gene_id)
tpm_data <- tpm_data[keep_ids, ]

# ─── 4. Log₂(TPM+1) transformation ───────────────
expr_log <- log2(tpm_data + 1)

## --- after you have built the vector `group` -----------------
meta$Group <- factor(group, levels = c("Control", "CAD"))

## PCA data frame ---------------------------------------------------------
pca_obj <- prcomp(t(expr_log), scale. = TRUE)
pca_df <- as.data.frame(pca_obj$x[, 1:2])           # PC1 & PC2 scores
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Group <- meta$Group                          # add group label
pca_df$Group <- factor(pca_df$Group, levels = rev(levels(pca_df$Group)))


# coloured by biology
pca_plot = ggplot(pca_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 3) +
  labs(
       x = "PC1 (32%)",
       y = "PC2 (15%)") +
  theme_minimal()



set.seed(123)
ump <- umap::umap(t(expr_log))      # or t(expr_mat_bc)

umap_df <- data.frame(
  UMAP1 = ump$layout[, 1],
  UMAP2 = ump$layout[, 2],
  Group = meta$Group
)

umap_df$Group <- factor(umap_df$Group, levels = rev(levels(umap_df$Group)))


umap_plot = ggplot(umap_df, aes(UMAP1, UMAP2, colour = Group)) +
  geom_point(size = 3) +
  theme_minimal()


pca_plot + umap_plot + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect')

ggsave("Figures/Figure1.png")

# ─── 5. Linear modeling with limma ───────────────
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(CADvsControl = CAD - Control, levels = design)

# 3. Continue the pipeline on the corrected matrix
fit  <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

fit <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# ─── 6. Extract & filter DEGs ────────────────────
res <- topTable(fit2, coef = "CADvsControl", number = Inf,
                adjust.method = "fdr", sort.by = "B")
res$ensembl_gene_id <- rownames(res)

# Apply cutoffs
lfc_cut <- 1.0
fdr_cut <- 0.05
degs <- res %>% filter(abs(logFC) >= lfc_cut, adj.P.Val < fdr_cut)
dim(degs)

# Save outputs
write.table(res, "data/TopGenes_GSE208194.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(degs, "data/DEGs_GSE208194.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ─── 7. Annotate DEGs with biotype ───────────────
biotype_info <- getBM(
  attributes  = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = degs$ensembl_gene_id,
  mart       = ensembl
)

degs_annot <- left_join(degs, biotype_info, by = "ensembl_gene_id") %>%
  filter(external_gene_name != "")


dim(degs_annot)
table(degs_annot$gene_biotype)

write.table(degs_annot, "data/degs_annotated_gse208194.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# ─── 8. Separate lncRNAs vs mRNAs ────────────────
lnc_df <- degs_annot %>% 
  filter(gene_biotype %in% c("lncRNA","lincRNA","antisense","processed_transcript"))

mrna_df <- degs_annot %>% 
  filter(gene_biotype == "protein_coding")

lnc_df

mrna_df

# ─── 9. Pathway enrichment on mRNAs ──────────────
entrez <- bitr(mrna_df$external_gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
kegg_res <- enrichKEGG(gene         = entrez$ENTREZID,
                       organism     = 'hsa',
                       pAdjustMethod= 'BH',
                       pvalueCutoff = 0.05)

# KEGGS analysis
as.data.frame(kegg_res)

p1 = enrichplot::dotplot(kegg_res, showCategory=100) + 
  ggtitle("KEGG Enrichment")


# Perform GO enrichment analysis for Biological Process
go_bp <- enrichGO(gene = entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

# Perform GO enrichment analysis for Biological Process
as.data.frame(go_bp)

# Visualize Biological Process enrichment results
p2 = enrichplot::dotplot(go_bp)

# Perform GO enrichment analysis for Cellular Component
go_cc <- enrichGO(gene = entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

# Perform GO enrichment analysis for Cellular Component
as.data.frame(go_cc)

# Visualize Cellular Component enrichment results
p3 = enrichplot::dotplot(go_cc)

# Perform GO enrichment analysis for Molecular Function
go_mf <- enrichGO(gene = entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

# Perform GO enrichment analysis for Molecular Function
as.data.frame(go_mf)

# Visualize Cellular Component enrichment results
p4 = enrichplot::dotplot(go_mf)

# ─── 10. Volcano plot ────────────────────────────

(p1+p2) / (p3+p4)+ plot_annotation(tag_levels = 'A') 
ggsave("Figures/Figure3.png", width = 15, height =12)


adjPVal_cutoff = 0.05
logFC_cutoff = 1

down <- rownames(res)[res$adj.P.Val < adjPVal_cutoff & res$logFC < -logFC_cutoff]
up <- rownames(res)[res$adj.P.Val < adjPVal_cutoff & res$logFC > logFC_cutoff]

## Print volcano plots with up genes in purple and down genes in green
keyvals.colour1 <- ifelse(
  rownames(res) %in% up, 'red',
  ifelse(rownames(res) %in% down, 'blue',
         'black'))

keyvals.colour1[is.na(keyvals.colour1)] <- 'black'
names(keyvals.colour1)[keyvals.colour1 == 'black'] <- 'Not significant'
names(keyvals.colour1)[keyvals.colour1 == 'red'] <- 'Up-regulated'
names(keyvals.colour1)[keyvals.colour1 == 'blue'] <- 'Down-regulated'


volcano_plot <-  EnhancedVolcano(res,
                       lab = NA,
                       x = 'logFC',
                       y = 'adj.P.Val',
                       title = NULL,  
                       subtitle = NULL, 
                       caption = NULL,
                       pCutoff = adjPVal_cutoff,
                       FCcutoff = logFC_cutoff,
                       pointSize = c(ifelse(res$adj.P.Val < adjPVal_cutoff & abs(res$logFC) > logFC_cutoff, 4, 2)),
                       colCustom = keyvals.colour1,
                       legendLabSize = 12,
                       legendIconSize = 8,
                       axisLabSize = 14,
                       legendPosition = 'top')

ggsave("Figures/Figure2.png")


