# =====================================================================
# Elastic Net Diagnostic Model with Probe Filtering,
# and Near-Zero Variance Removal
# =====================================================================

# 0. Set working directory
setwd("~/GitHub/cad")

# 1. Load libraries -----------------------------------------------------
library(data.table)    # fast data I/O
library(dplyr)         # data wrangling
library(ggplot2)       # plotting
library(caret)         # modeling framework
library(glmnet)        # elastic net / lasso / ridge
library(pROC)          # ROC/AUC calculations

# 2. Read feature list and expression data -----------------------------
best_probes   <- fread("data/best_probes.txt", stringsAsFactors = FALSE) %>% as.data.frame()
gset_exp_disc <- fread("data/gse208194_disc_expr.txt") %>% as.data.frame()
gset_exp_val  <- fread("data/gse113079_val_expr.txt")  %>% as.data.frame()

# 3. Filter probes present in both discovery & validation -------------
common_disc <- intersect(best_probes$ensembl_gene_id, colnames(gset_exp_disc))
common_val  <- intersect(best_probes$PROBE_ID,        colnames(gset_exp_val))
if (nrow(best_probes) > length(common_disc) || nrow(best_probes) > length(common_val)) {
  warning("Some probes dropped: not found in both datasets")
}
best_probes <- best_probes %>%
  filter(ensembl_gene_id %in% common_disc,
         PROBE_ID        %in% common_val)

# 4. Prepare discovery dataset ------------------------------------------
discovery <- gset_exp_disc[, c(best_probes$ensembl_gene_id, "group")]
colnames(discovery) <- c(best_probes$GeneSymbol, "group")
discovery$group <- factor(discovery$group, levels = c("Control","CAD"))

# 5. Prepare validation dataset -----------------------------------------
validation <- gset_exp_val[, c(best_probes$PROBE_ID, "group")]
colnames(validation) <- c(best_probes$GeneSymbol, "group")
validation$group <- factor(validation$group, levels = c("Control","CAD"))

# 6. Remove near-zero variance predictors from discovery & validation ---
nzv <- nearZeroVar(discovery %>% dplyr::select(-group))
if (length(nzv) > 0) {
  drop_genes <- names(discovery)[nzv]
  discovery  <- discovery  %>% dplyr::select(-all_of(drop_genes))
  validation <- validation %>% dplyr::select(-all_of(drop_genes))
  best_probes <- best_probes %>%
    filter(!GeneSymbol %in% drop_genes)
  message("Dropped near-zero variance genes: ", paste(drop_genes, collapse = ", "))
}

# 7. Split predictors (X) and outcome (y) --------------------------------
x_disc <- discovery %>% dplyr::select(-group)
y_disc <- discovery$group
x_val  <- validation  %>% dplyr::select(-group)
y_val  <- validation$group

# 8. Set seed for reproducibility ----------------------------------------
set.seed(42)

# 9. Define trainControl -----------------------------
ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 10,
  repeats         = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "all",
  returnResamp    = "all" 
)

# 10. Define tuning grid -------------------------------------------------
alpha_vals  <- seq(0, 1, by = 0.05)
lambda_vals <- 10^seq(-4, 2, length = 100) #lambda_vals <- 10^seq(-6, -1, length = 100)
tune_grid   <- expand.grid(alpha = alpha_vals, lambda = lambda_vals)

# 11. Train Elastic Net --------------------------------------------------
fit_enet <- train(
  x          = x_disc,
  y          = y_disc,
  method     = "glmnet",
  metric     = "ROC",
  trControl  = ctrl,
  tuneGrid   = tune_grid,
  preProcess = c("center","scale")
)

# 12. Review best tuning parameters and coefficients --------------------
cat("\nBest tuning parameters:\n")
print(fit_enet$bestTune)

coef_final <- coef(fit_enet$finalModel, s = fit_enet$bestTune$lambda)
coef_mat   <- as.matrix(coef_final)
cat("\nNon-zero coefficients:\n")
print(coef_mat[coef_mat[,1] != 0, , drop = FALSE])

plot(fit_enet)

## Coefficient Plot ###

library(ggplot2)
library(dplyr)
library(tidyr)

# Extract coefficients across lambdas
coef_path <- as.matrix(fit_enet$finalModel$beta)  # Features x Lambda
lambdas   <- fit_enet$finalModel$lambda           # Vector of lambdas

# Correct: transpose so that each row is a lambda
coef_df <- as.data.frame(t(coef_path))             # Lambda x Features
coef_df$lambda <- lambdas                         # Add lambda as a column

# Reshape to long format
coef_long <- coef_df %>%
  pivot_longer(
    cols = -lambda,
    names_to = "feature",
    values_to = "coefficient"
  )

# Plot
p2  = ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = feature)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Log(Lambda)", y = "Coefficient") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
        )  # optional: hide legend if too crowded


### Lambda Plot ####

library(ggplot2)
library(dplyr)

# 1. Extract prediction data
pred_df <- fit_enet$pred

# 2. Focus only on best alpha
best_alpha <- fit_enet$bestTune$alpha
pred_df_best_alpha <- pred_df %>%
  filter(alpha == best_alpha)

# 3. Find positive class name
positive_class <- levels(pred_df$obs)[2]

# 4. Calculate MSE and Standard Error
mse_by_lambda <- pred_df_best_alpha %>%
  group_by(lambda) %>%
  summarize(
    mse = mean((ifelse(obs == positive_class, 1, 0) - .data[[positive_class]])^2),
    se  = sd((ifelse(obs == positive_class, 1, 0) - .data[[positive_class]])^2) / sqrt(n()),
    n = n()
  ) %>%
  ungroup()

# 5. Find best lambda (smallest mean MSE)
best_lambda_mse <- mse_by_lambda %>%
  filter(mse == min(mse)) %>%
  slice(1) %>%
  pull(lambda)

# 6. Plot with error bars
p1 <- ggplot(mse_by_lambda, aes(x = log(lambda), y = mse)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = mse - se, ymax = mse + se), width = 0.05) +
  geom_vline(xintercept = log(best_lambda_mse), linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = log(best_lambda_mse), y = max(mse_by_lambda$mse), 
           label = "Best Lambda", hjust = 1.5, vjust = 1.5, size = 3.5, color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Log(Lambda)",
    y = "Mean Squared Error (MSE)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


library(patchwork)
p1+p2+ plot_annotation(tag_levels = 'A') 

ggsave("Figures/Figure4.png")



# =====================================================================
# 12 b.  Bootstrap feature-selection stability (1000 resamples)
# =====================================================================

set.seed(123)                       # reproducibility
B <- 1000                            # number of bootstraps
alpha_best   <- fit_enet$bestTune$alpha
lambda_best  <- fit_enet$bestTune$lambda

# matrix that will store 1 = "selected" for every gene in every bootstrap
sel_mat <- matrix(0, nrow = B, ncol = ncol(x_disc))
colnames(sel_mat) <- colnames(x_disc)

# convert outcome to 0/1 once so we can reuse it
y_num <- ifelse(y_disc == "CAD", 1, 0)

for (b in seq_len(B)) {
  idx   <- sample.int(nrow(x_disc), replace = TRUE)   # bootstrap rows
  x_boot <- as.matrix(x_disc[idx, ])
  y_boot <- y_num[idx]
  
  fit_b <- glmnet(
    x      = x_boot,
    y      = y_boot,
    alpha  = alpha_best,
    lambda = lambda_best,
    family = "binomial",
    standardize = TRUE
  )
  
  coef_b <- coef(fit_b)[, 1]                    # vector of coefficients
  sel_genes <- names(coef_b[coef_b != 0])
  sel_genes <- setdiff(sel_genes, "(Intercept)")
  if (length(sel_genes))
    sel_mat[b, sel_genes] <- 1
}

# ---------------------------------------------------------------------
# Selection frequency table
freq <- colMeans(sel_mat)           # proportion of times selected
freq_df <- data.frame(
  Gene = names(freq),
  SelectionFreq = round(100 * freq, 1)          # %
) %>% arrange(desc(SelectionFreq))

# Genes selected in ≥70 % of bootstraps
stable_genes <- subset(freq_df, SelectionFreq >= 50)

cat("\nElastic-net feature-selection stability (≥50 % of 200 bootstraps):\n")
print(stable_genes, row.names = FALSE)

# Save results
write.table(freq_df, "data/ENet_bootstrap_selection_frequency.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

## Optional: bar plot of the top 20 most stable genes -------------------
library(ggplot2)
top20 <- head(freq_df, 20)

ggplot(top20, aes(x = reorder(Gene, SelectionFreq), y = SelectionFreq)) +
  geom_col() +
  coord_flip() +
  labs(
       x = NULL, y = "Selection frequency (%)") +
  theme_minimal()

ggsave("Figures/Figure6.png", width = 6, height = 4)


# 13. Predict probabilities on validation set ----------------------------
pred_prob <- predict(fit_enet, newdata = x_val, type = "prob")

# 14. Compute ROC & AUC --------------------------------------------------
library(pROC)

# Compute ROC and AUC
roc_obj <- roc(
  response  = y_val,
  predictor = pred_prob$CAD,
  levels    = c("Control", "CAD")
)
auc_val <- auc(roc_obj)

# Compute 95% CI for AUC
ci_auc <- ci.auc(roc_obj)

png("Figures/Figure5.png", width = 800, height = 600, res = 120)  # Open a PNG device


# Plot ROC curve
plot(
  roc_obj,
  col = "#1c61b6",
  legacy.axes = TRUE

)

# Add AUC and CI text
text(
  x = 0.55, y = 0.1,  # Lower right position
  labels = paste0("AUC = ", round(auc_val, 3), 
                  " (95% CI: ", round(ci_auc[1], 3), "-", round(ci_auc[3], 3), ")"),
  adj = 0
)

dev.off()  # Close the PNG device


# 15. Determine optimal decision threshold -------------------------------
opt <- coords(
  roc_obj, "best",
  ret         = c("threshold","sensitivity","specificity"),
  best.method = "youden"
)
cat("\nOptimal threshold:", round(opt[["threshold"]], 3),
    "Sensitivity:", round(opt[["sensitivity"]], 3),
    "Specificity:", round(opt[["specificity"]], 3), "\n")

pred_class_opt <- factor(
  ifelse(pred_prob$CAD > opt[["threshold"]], "CAD", "Control"),
  levels = levels(y_val)
)
cm_opt <- confusionMatrix(pred_class_opt, y_val, positive = "CAD")
cat("\nConfusion Matrix at optimal threshold:\n")
print(cm_opt)

# 16. Probability calibration (Platt scaling) -----------------
cal_mod <- glm(as.numeric(y_val == "CAD") ~ pred_prob$CAD, family = binomial)
pred_prob_cal <- predict(
  cal_mod,
  newdata = data.frame(pred_prob.CAD = pred_prob$CAD),
  type    = "response"
)
roc_cal <- roc(
  response  = y_val,
  predictor = pred_prob_cal,
  levels    = c("Control","CAD")
)
auc_cal <- auc(roc_cal)
cat("\nCalibrated AUC:", round(auc_cal, 3), "\n")
plot(
  roc_cal,
  main       = sprintf("Calibrated ROC (AUC = %.3f)", auc_cal),
  col        = "#d95f02",
  legacy.axes= TRUE
)

opt_cal <- coords(
  roc_cal, "best",
  ret         = c("threshold","sensitivity","specificity"),
  best.method = "youden"
)
pred_class_cal <- factor(
  ifelse(pred_prob_cal > opt_cal[["threshold"]], "CAD", "Control"),
  levels = levels(y_val)
)
cm_cal <- confusionMatrix(pred_class_cal, y_val, positive = "CAD")
cat("\nConfusion Matrix after calibration:\n")
print(cm_cal)


best      <- fit_enet$bestTune
idx       <- with(fit_enet$resample,
                  alpha == best$alpha & lambda == best$lambda)
cv_auc    <- mean(fit_enet$resample$ROC[idx])   # ≈ 0.8396
cv_auc

