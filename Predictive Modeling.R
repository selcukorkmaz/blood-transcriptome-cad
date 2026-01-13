# =====================================================================
# Elastic Net Diagnostic Model with Guarded Resampling (fastml)
# Leakage-aware DEG selection and nested CV
# =====================================================================

# 0. Set working directory
setwd("~/GitHub/cad")

# 1. Load libraries -----------------------------------------------------
library(data.table)    # fast data I/O
library(dplyr)         # data wrangling
library(ggplot2)       # plotting
library(fastml)        # guarded resampling + tuning
library(rsample)       # resampling splits
library(recipes)       # preprocessing
library(glmnet)        # elastic net / lasso / ridge
library(limma)         # DEG selection within folds
library(pROC)          # ROC/AUC calculations
library(tidyr)         # pivot_longer
library(yardstick)     # confusion matrix
library(tune)          # extract_fit_parsnip
library(patchwork)     # plot composition

# 2. Read feature list and expression data -----------------------------
best_probes   <- fread("data/best_probes.txt", stringsAsFactors = FALSE) %>% as.data.frame()
gset_exp_disc <- fread("data/gse208194_disc_expr.txt") %>% as.data.frame()
gset_exp_val  <- fread("data/gse113079_val_expr.txt")  %>% as.data.frame()

# 3. Filter probes present in both discovery & validation --------------
common_disc <- intersect(best_probes$ensembl_gene_id, colnames(gset_exp_disc))
common_val  <- intersect(best_probes$PROBE_ID,        colnames(gset_exp_val))
if (nrow(best_probes) > length(common_disc) || nrow(best_probes) > length(common_val)) {
  warning("Some probes dropped: not found in both datasets")
}
best_probes <- best_probes %>%
  filter(ensembl_gene_id %in% common_disc,
         PROBE_ID        %in% common_val)

# 4. Prepare discovery dataset (global feature mapping) ----------------
discovery <- gset_exp_disc[, c(best_probes$ensembl_gene_id, "group")]
colnames(discovery) <- c(best_probes$GeneSymbol, "group")
discovery$group <- factor(discovery$group, levels = c("Control","CAD"))

# 5. Prepare validation dataset ----------------------------------------
validation <- gset_exp_val[, c(best_probes$PROBE_ID, "group")]
colnames(validation) <- c(best_probes$GeneSymbol, "group")
validation$group <- factor(validation$group, levels = c("Control","CAD"))

# 6. Full discovery matrix for leakage-free CV -------------------------
# NOTE: This is intentionally NOT filtered by global DEGs.
discovery_full <- gset_exp_disc
discovery_full$group <- factor(discovery_full$group, levels = c("Control","CAD"))
numeric_cols <- names(discovery_full)[sapply(discovery_full, is.numeric)]
numeric_cols <- setdiff(numeric_cols, "group")
discovery_full <- discovery_full[, c(numeric_cols, "group"), drop = FALSE]

# 7. Modeling parameters ------------------------------------------------
set.seed(42)

lfc_cut <- 1.0
fdr_cut <- 0.05
min_degs <- 10

outer_folds   <- 10
inner_folds   <- 10
inner_repeats <- 5

alpha_vals  <- seq(0, 1, by = 0.05)
lambda_vals_log10 <- seq(-4, 2, length = 100)
lambda_vals <- 10^lambda_vals_log10

tune_grid <- list(
  logistic_reg = list(
    glmnet = list(
      penalty = lambda_vals_log10,
      mixture = alpha_vals
    )
  )
)

# 8. Helper functions ---------------------------------------------------
compute_degs <- function(train_df, lfc_cut, fdr_cut, min_degs = 10) {
  x_train <- train_df %>% dplyr::select(-group)
  y_train <- train_df$group

  if (ncol(x_train) == 0 || length(unique(y_train)) < 2) {
    return(character(0))
  }

  expr_mat <- t(as.matrix(x_train))
  design <- model.matrix(~0 + y_train)
  colnames(design) <- levels(y_train)
  contrast <- makeContrasts(CADvsControl = CAD - Control, levels = design)

  fit <- lmFit(expr_mat, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = "CADvsControl", number = Inf,
                  adjust.method = "fdr", sort.by = "B")
  res$gene_id <- rownames(res)

  degs <- res %>% dplyr::filter(abs(logFC) >= lfc_cut, adj.P.Val < fdr_cut)
  if (nrow(degs) < min_degs) {
    degs <- res %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = min(min_degs, nrow(res)))
  }
  degs$gene_id
}

extract_prob <- function(pred_df, positive_class) {
  if (is.list(pred_df) && !inherits(pred_df, "data.frame")) {
    pred_df <- pred_df[[1]]
  }
  if (!inherits(pred_df, "data.frame")) {
    stop("Unexpected prediction format from fastml.")
  }

  prob_col <- paste0(".pred_", positive_class)
  if (!prob_col %in% names(pred_df)) {
    prob_col <- paste0(".pred_p", positive_class)
  }
  if (!prob_col %in% names(pred_df)) {
    prob_cols <- grep("^\\.pred_", names(pred_df), value = TRUE)
    if (length(prob_cols) < 1) {
      stop("No probability column found in predictions.")
    }
    prob_col <- prob_cols[1]
  }

  list(truth = pred_df$truth, prob = pred_df[[prob_col]])
}

fit_fastml_fold <- function(train_df, test_df, tune_grid, inner_folds, inner_repeats, seed) {
  rec <- recipe(group ~ ., data = train_df) %>%
    step_nzv(all_predictors()) %>%
    step_normalize(all_predictors())

  fastml(
    train_data = train_df,
    test_data = test_df,
    label = "group",
    algorithms = "logistic_reg",
    algorithm_engines = list(logistic_reg = "glmnet"),
    resampling_method = "repeatedcv",
    folds = inner_folds,
    repeats = inner_repeats,
    tune_params = tune_grid,
    use_default_tuning = FALSE,
    tuning_strategy = "grid",
    metric = "roc_auc",
    recipe = rec,
    event_class = "second",
    seed = seed,
    bootstrap_ci = FALSE,
    verbose = FALSE
  )
}

run_cv_pipeline <- function(data, outer_resamples, deg_strategy, lfc_cut, fdr_cut, min_degs,
                            tune_grid, inner_folds, inner_repeats, seed) {
  positive_class <- levels(data$group)[2]
  deg_sets <- vector("list", length(outer_resamples$splits))
  fold_auc <- rep(NA_real_, length(outer_resamples$splits))
  all_truth <- character(0)
  all_prob <- numeric(0)

  global_degs <- NULL
  if (deg_strategy == "global") {
    global_degs <- compute_degs(data, lfc_cut, fdr_cut, min_degs)
  }

  for (i in seq_along(outer_resamples$splits)) {
    split <- outer_resamples$splits[[i]]
    train_df <- rsample::analysis(split)
    test_df  <- rsample::assessment(split)

    if (length(unique(train_df$group)) < 2 || length(unique(test_df$group)) < 2) {
      warning("Skipped fold ", i, ": class imbalance in split.")
      next
    }

    degs <- if (deg_strategy == "global") {
      global_degs
    } else {
      compute_degs(train_df, lfc_cut, fdr_cut, min_degs)
    }

    degs <- setdiff(intersect(degs, colnames(train_df)), "group")
    deg_sets[[i]] <- degs

    if (length(degs) < 2) {
      warning("Skipped fold ", i, ": too few DEGs after filtering.")
      next
    }

    train_sel <- train_df[, c(degs, "group"), drop = FALSE]
    test_sel  <- test_df[,  c(degs, "group"), drop = FALSE]

    set.seed(seed + i)
    fit_fold <- fit_fastml_fold(train_sel, test_sel, tune_grid, inner_folds, inner_repeats, seed + i)
    pred_info <- extract_prob(fit_fold$predictions, positive_class)

    truth <- factor(pred_info$truth, levels = levels(data$group))
    prob  <- as.numeric(pred_info$prob)

    all_truth <- c(all_truth, as.character(truth))
    all_prob  <- c(all_prob, prob)

    roc_fold <- pROC::roc(response = truth, predictor = prob, levels = levels(data$group))
    fold_auc[i] <- as.numeric(pROC::auc(roc_fold))
  }

  all_truth <- factor(all_truth, levels = levels(data$group))
  roc_all <- pROC::roc(response = all_truth, predictor = all_prob, levels = levels(data$group))
  auc_all <- as.numeric(pROC::auc(roc_all))
  ci_all <- pROC::ci.auc(roc_all)

  list(
    auc = auc_all,
    ci = ci_all,
    fold_auc = fold_auc,
    deg_sets = deg_sets
  )
}

# 9. Leakage audit: baseline vs corrected CV ---------------------------
outer_resamples <- rsample::vfold_cv(discovery_full, v = outer_folds, strata = group)

discovery_full[1:5,1:5]
dim(discovery_full)

cv_clean <- run_cv_pipeline(
  data = discovery_full,
  outer_resamples = outer_resamples,
  deg_strategy = "per_fold",
  lfc_cut = lfc_cut,
  fdr_cut = fdr_cut,
  min_degs = min_degs,
  tune_grid = tune_grid,
  inner_folds = inner_folds,
  inner_repeats = inner_repeats,
  seed = 42
)


cat("\nLeakage-corrected CV AUC (fold-wise DEG selection): ",
    round(cv_clean$auc, 3),
    " (95% CI: ", round(cv_clean$ci[1], 3), "-",
    round(cv_clean$ci[3], 3), ")\n", sep = "")

# Optional diagnostics: DEG counts and stability
deg_counts <- sapply(cv_clean$deg_sets, length)
deg_count_df <- data.frame(
  Fold = seq_along(deg_counts),
  DEG_Count = deg_counts
)
write.table(deg_count_df, "data/DEG_counts_per_fold.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

deg_freq <- sort(table(unlist(cv_clean$deg_sets)) / length(cv_clean$deg_sets), decreasing = TRUE)
deg_freq_df <- data.frame(
  Gene = names(deg_freq),
  SelectionFreq = round(100 * as.numeric(deg_freq), 1)
)
write.table(deg_freq_df, "data/DEG_selection_frequency_cv.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# 10. Train final model on full discovery and validate ------------------
rec_final <- recipe(group ~ ., data = discovery) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors())

fit_enet <- fastml(
  train_data = discovery,
  test_data = validation,
  label = "group",
  algorithms = "logistic_reg",
  algorithm_engines = list(logistic_reg = "glmnet"),
  resampling_method = "repeatedcv",
  folds = inner_folds,
  repeats = inner_repeats,
  tune_params = tune_grid,
  use_default_tuning = FALSE,
  tuning_strategy = "grid",
  metric = "roc_auc",
  recipe = rec_final,
  event_class = "second",
  seed = 42,
  bootstrap_ci = FALSE,
  verbose = FALSE
)

get_param_value <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }
  if (inherits(x, "quosure")) {
    x <- rlang::eval_tidy(x)
  }
  as.numeric(x)[1]
}

best_workflow <- fit_enet$best_model[[1]]
best_parsnip <- tune::extract_fit_parsnip(best_workflow)
best_params <- best_parsnip$spec$args

alpha_best <- get_param_value(best_params$mixture, default = 1)
lambda_best <- get_param_value(best_params$penalty, default = NA_real_)

cat("\nBest tuning parameters (fastml):\n")
print(list(alpha = alpha_best, lambda = lambda_best))

# Preprocess training data for plotting and stability analyses
train_processed <- bake(fit_enet$preprocessor, new_data = discovery)
x_disc <- as.matrix(train_processed %>% dplyr::select(-group))
y_num <- ifelse(train_processed$group == "CAD", 1, 0)

cv_fit <- cv.glmnet(
  x_disc,
  y_num,
  family = "binomial",
  alpha = alpha_best,
  lambda = lambda_vals,
  standardize = FALSE
)

if (!is.finite(lambda_best)) {
  lambda_best <- cv_fit$lambda.min
}

glmnet_path <- glmnet(
  x_disc,
  y_num,
  family = "binomial",
  alpha = alpha_best,
  lambda = lambda_vals,
  standardize = FALSE
)

coef_final <- coef(glmnet_path, s = lambda_best)
coef_mat <- as.matrix(coef_final)

cat("\nNon-zero coefficients:\n")
print(coef_mat[coef_mat[, 1] != 0, , drop = FALSE])

# Coefficient Plot ------------------------------------------------------
coef_path <- as.matrix(glmnet_path$beta)
lambdas   <- glmnet_path$lambda

coef_df <- as.data.frame(t(coef_path))
coef_df$lambda <- lambdas

coef_long <- coef_df %>%
  pivot_longer(
    cols = -lambda,
    names_to = "feature",
    values_to = "coefficient"
  )

p2 <- ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = feature)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Log(Lambda)", y = "Coefficient") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Lambda Plot -----------------------------------------------------------
cv_df <- data.frame(
  lambda = cv_fit$lambda,
  cvm = cv_fit$cvm,
  cvsd = cv_fit$cvsd
)

best_lambda_mse <- cv_fit$lambda.min

p1 <- ggplot(cv_df, aes(x = log(lambda), y = cvm)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), width = 0.05) +
  geom_vline(xintercept = log(best_lambda_mse), linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = log(best_lambda_mse), y = max(cv_df$cvm),
           label = "Best Lambda", hjust = 1.5, vjust = 1.5, size = 3.5, color = "red") +
  theme_minimal(base_size = 14) +
  labs(x = "Log(Lambda)", y = "Mean CV Deviance") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

p1 + p2 + plot_annotation(tag_levels = "A")
ggsave("Figures/Figure4.png")

# =====================================================================
# 12 b.  Bootstrap feature-selection stability (1000 resamples)
# =====================================================================

set.seed(123)
B <- 1000

sel_mat <- matrix(0, nrow = B, ncol = ncol(x_disc))
colnames(sel_mat) <- colnames(x_disc)

for (b in seq_len(B)) {
  idx <- sample.int(nrow(x_disc), replace = TRUE)
  x_boot <- x_disc[idx, , drop = FALSE]
  y_boot <- y_num[idx]

  fit_b <- glmnet(
    x = x_boot,
    y = y_boot,
    alpha = alpha_best,
    lambda = lambda_best,
    family = "binomial",
    standardize = FALSE
  )

  coef_b <- coef(fit_b)[, 1]
  sel_genes <- names(coef_b[coef_b != 0])
  sel_genes <- setdiff(sel_genes, "(Intercept)")
  if (length(sel_genes)) {
    sel_mat[b, sel_genes] <- 1
  }
}

freq <- colMeans(sel_mat)
freq_df <- data.frame(
  Gene = names(freq),
  SelectionFreq = round(100 * freq, 1)
) %>% arrange(desc(SelectionFreq))

stable_genes <- subset(freq_df, SelectionFreq >= 50)

cat("\nElastic-net feature-selection stability (>= 50% of bootstraps):\n")
print(stable_genes, row.names = FALSE)

write.table(freq_df, "data/ENet_bootstrap_selection_frequency.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

top20 <- head(freq_df, 20)

ggplot(top20, aes(x = reorder(Gene, SelectionFreq), y = SelectionFreq)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Selection frequency (%)") +
  theme_minimal()

ggsave("Figures/Figure6.png", width = 6, height = 4)

# 13. Predict probabilities on validation set --------------------------
pred_df <- fit_enet$predictions
positive_class <- levels(validation$group)[2]
pred_info <- extract_prob(pred_df, positive_class)
pred_prob <- pred_info$prob
y_val <- validation$group

# 14. Compute ROC & AUC -------------------------------------------------
roc_obj <- roc(
  response = y_val,
  predictor = pred_prob,
  levels = c("Control", "CAD")
)
auc_val <- auc(roc_obj)
ci_auc <- ci.auc(roc_obj)

png("Figures/Figure5.png", width = 800, height = 600, res = 120)

plot(
  roc_obj,
  col = "#1c61b6",
  legacy.axes = TRUE
)

text(
  x = 0.55, y = 0.1,
  labels = paste0("AUC = ", round(auc_val, 3),
                  " (95% CI: ", round(ci_auc[1], 3), "-", round(ci_auc[3], 3), ")"),
  adj = 0
)

dev.off()

# 15. Determine optimal decision threshold ------------------------------
opt <- coords(
  roc_obj, "best",
  ret = c("threshold", "sensitivity", "specificity"),
  best.method = "youden"
)
cat("\nOptimal threshold:", round(opt[["threshold"]], 3),
    "Sensitivity:", round(opt[["sensitivity"]], 3),
    "Specificity:", round(opt[["specificity"]], 3), "\n")

pred_class_opt <- factor(
  ifelse(pred_prob > opt[["threshold"]], "CAD", "Control"),
  levels = levels(y_val)
)

cm_opt <- yardstick::conf_mat(
  data.frame(truth = y_val, estimate = pred_class_opt),
  truth = truth,
  estimate = estimate
)

cat("\nConfusion Matrix at optimal threshold:\n")
print(cm_opt)

# 16. Probability calibration (Platt scaling) ---------------------------
cal_mod <- glm(as.numeric(y_val == "CAD") ~ pred_prob, family = binomial)
pred_prob_cal <- predict(
  cal_mod,
  newdata = data.frame(pred_prob = pred_prob),
  type = "response"
)

roc_cal <- roc(
  response = y_val,
  predictor = pred_prob_cal,
  levels = c("Control", "CAD")
)
auc_cal <- auc(roc_cal)
cat("\nCalibrated AUC:", round(auc_cal, 3), "\n")

plot(
  roc_cal,
  main = sprintf("Calibrated ROC (AUC = %.3f)", auc_cal),
  col = "#d95f02",
  legacy.axes = TRUE
)

opt_cal <- coords(
  roc_cal, "best",
  ret = c("threshold", "sensitivity", "specificity"),
  best.method = "youden"
)
pred_class_cal <- factor(
  ifelse(pred_prob_cal > opt_cal[["threshold"]], "CAD", "Control"),
  levels = levels(y_val)
)

cm_cal <- yardstick::conf_mat(
  data.frame(truth = y_val, estimate = pred_class_cal),
  truth = truth,
  estimate = estimate
)

cat("\nConfusion Matrix after calibration:\n")
print(cm_cal)
