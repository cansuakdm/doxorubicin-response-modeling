############################################################
# RANDOM FOREST RFE WITH 5-FOLD CV FOR ALL GENE SETS
############################################################

setwd("~/Downloads/DepMap_Project") 

library(data.table)
library(randomForest)

final_dataset <- fread("Dataset_A_gene_expression_ic50.csv")
X <- fread("filtered_gene_matrix.csv")

y <- final_dataset$log_ic50
X <- as.data.table(X)

cat("Number of samples =", nrow(X), "\n")
cat("Number of selected genes =", ncol(X), "\n")

calculate_metrics <- function(actual, predicted) {
  mae <- mean(abs(predicted - actual))
  rmse <- sqrt(mean((predicted - actual)^2))
  r2 <- cor(predicted, actual)^2
  
  return(c(MAE = mae, RMSE = rmse, R2 = r2))
}

############################################################
# PART 1: CREATE RFE GENE SETS
############################################################

set.seed(42)

train_idx <- sample(
  1:nrow(X),
  size = 0.8 * nrow(X)
)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

y_train <- y[train_idx]
y_test  <- y[-train_idx]

current_genes <- colnames(X_train)

rfe_results <- data.frame(
  Stage = character(),
  iteration = integer(),
  n_genes = integer(),
  MAE = numeric(),
  RMSE = numeric(),
  R2 = numeric()
)

gene_list <- list()

############################################################
# PART 2: RFE BY REMOVING LOWEST 3% IMPORTANT GENES
############################################################

iteration <- 1

while (length(current_genes) > 52) {
  
  cat(
    "3% RFE iteration:",
    iteration,
    "| Number of genes:",
    length(current_genes),
    "\n"
  )
  
  set.seed(42)
  
  rf_model <- randomForest(
    x = X_train[, current_genes, with = FALSE],
    y = y_train,
    ntree = 500,
    importance = TRUE
  )
  
  pred <- predict(
    rf_model,
    X_test[, current_genes, with = FALSE]
  )
  
  metrics <- calculate_metrics(
    actual = y_test,
    predicted = pred
  )
  
  rfe_results <- rbind(
    rfe_results,
    data.frame(
      Stage = "RFE_3_percent",
      iteration = iteration,
      n_genes = length(current_genes),
      MAE = as.numeric(metrics["MAE"]),
      RMSE = as.numeric(metrics["RMSE"]),
      R2 = as.numeric(metrics["R2"])
    )
  )
  
  gene_list[[paste0("RFE_3_percent_", iteration)]] <- current_genes
  
  importance_vals <- importance(rf_model)
  
  importance_df <- data.frame(
    Gene = rownames(importance_vals),
    IncMSE = importance_vals[, "%IncMSE"]
  )
  
  importance_df <- importance_df[
    order(importance_df$IncMSE),
  ]
  
  remove_n <- ceiling(
    0.03 * nrow(importance_df)
  )
  
  genes_to_remove <- importance_df$Gene[1:remove_n]
  
  current_genes <- setdiff(
    current_genes,
    genes_to_remove
  )
  
  iteration <- iteration + 1
}

############################################################
# PART 3: REMOVE REMAINING GENES ONE BY ONE
############################################################

single_iteration <- 1

while (length(current_genes) > 1) {
  
  cat(
    "Single-gene elimination iteration:",
    single_iteration,
    "| Number of genes:",
    length(current_genes),
    "\n"
  )
  
  set.seed(42)
  
  rf_model <- randomForest(
    x = X_train[, current_genes, with = FALSE],
    y = y_train,
    ntree = 500,
    importance = TRUE
  )
  
  pred <- predict(
    rf_model,
    X_test[, current_genes, with = FALSE]
  )
  
  metrics <- calculate_metrics(
    actual = y_test,
    predicted = pred
  )
  
  rfe_results <- rbind(
    rfe_results,
    data.frame(
      Stage = "Single_gene_elimination",
      iteration = single_iteration,
      n_genes = length(current_genes),
      MAE = as.numeric(metrics["MAE"]),
      RMSE = as.numeric(metrics["RMSE"]),
      R2 = as.numeric(metrics["R2"])
    )
  )
  
  gene_list[[paste0("Single_gene_elimination_", single_iteration)]] <- current_genes
  
  importance_vals <- importance(rf_model)
  
  importance_df <- data.frame(
    Gene = rownames(importance_vals),
    IncMSE = importance_vals[, "%IncMSE"]
  )
  
  importance_df <- importance_df[
    order(importance_df$IncMSE),
  ]
  
  gene_to_remove <- importance_df$Gene[1]
  
  current_genes <- setdiff(
    current_genes,
    gene_to_remove
  )
  
  single_iteration <- single_iteration + 1
}

############################################################
# PART 4: 5-FOLD CV FUNCTION FOR ONE GENE SET
############################################################

run_rf_cv_for_gene_set <- function(X_data, y_data, gene_set, model_name) {
  
  set.seed(42)
  
  k <- 5
  
  folds <- sample(
    rep(1:k, length.out = nrow(X_data))
  )
  
  cv_results <- data.frame(
    Model = character(),
    Fold = integer(),
    n_genes = integer(),
    MAE = numeric(),
    RMSE = numeric(),
    R2 = numeric()
  )
  
  cv_predictions <- data.frame(
    Model = character(),
    Actual = numeric(),
    Predicted = numeric(),
    Fold = integer()
  )
  
  for (i in 1:k) {
    
    cat("Running", model_name, "- Fold", i, "\n")
    
    test_idx <- which(folds == i)
    train_idx <- setdiff(1:nrow(X_data), test_idx)
    
    X_train_cv <- X_data[train_idx, gene_set, with = FALSE]
    X_test_cv  <- X_data[test_idx, gene_set, with = FALSE]
    
    y_train_cv <- y_data[train_idx]
    y_test_cv  <- y_data[test_idx]
    
    set.seed(42)
    
    rf_cv_model <- randomForest(
      x = X_train_cv,
      y = y_train_cv,
      ntree = 500
    )
    
    pred_cv <- predict(
      rf_cv_model,
      X_test_cv
    )
    
    metrics_cv <- calculate_metrics(
      actual = y_test_cv,
      predicted = pred_cv
    )
    
    cv_results <- rbind(
      cv_results,
      data.frame(
        Model = model_name,
        Fold = i,
        n_genes = length(gene_set),
        MAE = as.numeric(metrics_cv["MAE"]),
        RMSE = as.numeric(metrics_cv["RMSE"]),
        R2 = as.numeric(metrics_cv["R2"])
      )
    )
    
    cv_predictions <- rbind(
      cv_predictions,
      data.frame(
        Model = model_name,
        Actual = y_test_cv,
        Predicted = pred_cv,
        Fold = i
      )
    )
  }
  
  overall_mae <- mean(
    abs(cv_predictions$Predicted - cv_predictions$Actual)
  )
  
  overall_rmse <- sqrt(
    mean((cv_predictions$Predicted - cv_predictions$Actual)^2)
  )
  
  overall_r2 <- cor(
    cv_predictions$Predicted,
    cv_predictions$Actual
  )^2
  
  cv_summary <- data.frame(
    Model = model_name,
    n_genes = length(gene_set),
    Mean_MAE = mean(cv_results$MAE),
    Mean_RMSE = mean(cv_results$RMSE),
    Mean_R2 = mean(cv_results$R2),
    Overall_MAE = overall_mae,
    Overall_RMSE = overall_rmse,
    Overall_R2 = overall_r2
  )
  
  return(
    list(
      fold_results = cv_results,
      predictions = cv_predictions,
      summary = cv_summary
    )
  )
}

############################################################
# PART 5: RUN 5-FOLD CV FOR ALL RFE GENE SETS
############################################################

all_cv_summary <- data.frame()
all_cv_folds <- data.frame()
all_cv_predictions <- data.frame()

for (gene_set_name in names(gene_list)) {
  
  gene_set <- gene_list[[gene_set_name]]
  
  cat(
    "Running 5-fold CV for",
    gene_set_name,
    "| Number of genes:",
    length(gene_set),
    "\n"
  )
  
  cv_output <- run_rf_cv_for_gene_set(
    X_data = X,
    y_data = y,
    gene_set = gene_set,
    model_name = gene_set_name
  )
  
  all_cv_summary <- rbind(
    all_cv_summary,
    cv_output$summary
  )
  
  all_cv_folds <- rbind(
    all_cv_folds,
    cv_output$fold_results
  )
  
  all_cv_predictions <- rbind(
    all_cv_predictions,
    cv_output$predictions
  )
}

############################################################
# PART 6: SELECT BEST MODELS
############################################################

best_rmse_row <- which.min(all_cv_summary$Overall_RMSE)
best_r2_row <- which.max(all_cv_summary$Overall_R2)

best_rmse_model <- all_cv_summary[best_rmse_row, ]
best_r2_model <- all_cv_summary[best_r2_row, ]

cat("\nBest model by Overall RMSE:\n")
print(best_rmse_model)

cat("\nBest model by Overall R2:\n")
print(best_r2_model)

############################################################
# PART 7: SAVE RESULTS
############################################################

write.csv(
  rfe_results,
  "rf_rfe_single_split_results.csv",
  row.names = FALSE
)

write.csv(
  all_cv_summary,
  "rf_rfe_all_gene_sets_5fold_cv_summary.csv",
  row.names = FALSE
)

write.csv(
  all_cv_folds,
  "rf_rfe_all_gene_sets_5fold_cv_folds.csv",
  row.names = FALSE
)

write.csv(
  all_cv_predictions,
  "rf_rfe_all_gene_sets_5fold_cv_predictions.csv",
  row.names = FALSE
)

best_rmse_genes <- gene_list[[best_rmse_model$Model]]
best_r2_genes <- gene_list[[best_r2_model$Model]]

write.csv(
  data.frame(gene = best_rmse_genes),
  "rf_best_genes_by_cv_rmse.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = best_r2_genes),
  "rf_best_genes_by_cv_r2.csv",
  row.names = FALSE
)

############################################################
# PART 8: PLOTS
############################################################

png(
  "rf_all_rfe_gene_sets_cv_RMSE.png",
  width = 2400,
  height = 1800,
  res = 300
)

plot(
  all_cv_summary$n_genes,
  all_cv_summary$Overall_RMSE,
  type = "b",
  pch = 19,
  xlab = "Number of Genes",
  ylab = "Overall CV RMSE",
  main = "Overall CV RMSE Across RFE Gene Sets"
)

points(
  all_cv_summary$n_genes[best_rmse_row],
  all_cv_summary$Overall_RMSE[best_rmse_row],
  pch = 19,
  cex = 2
)

text(
  all_cv_summary$n_genes[best_rmse_row],
  all_cv_summary$Overall_RMSE[best_rmse_row],
  labels = paste0(
    "Best: ",
    all_cv_summary$n_genes[best_rmse_row],
    " genes"
  ),
  pos = 3
)

dev.off()

png(
  "rf_all_rfe_gene_sets_cv_R2.png",
  width = 2400,
  height = 1800,
  res = 300
)

plot(
  all_cv_summary$n_genes,
  all_cv_summary$Overall_R2,
  type = "b",
  pch = 19,
  xlab = "Number of Genes",
  ylab = "Overall CV R²",
  main = "Overall CV R² Across RFE Gene Sets"
)

points(
  all_cv_summary$n_genes[best_r2_row],
  all_cv_summary$Overall_R2[best_r2_row],
  pch = 19,
  cex = 2
)

text(
  all_cv_summary$n_genes[best_r2_row],
  all_cv_summary$Overall_R2[best_r2_row],
  labels = paste0(
    "Best: ",
    all_cv_summary$n_genes[best_r2_row],
    " genes"
  ),
  pos = 3
)

dev.off()

############################################################