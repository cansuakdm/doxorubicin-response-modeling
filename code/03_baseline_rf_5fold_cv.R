############################################################
# BASELINE RANDOM FOREST 5-FOLD CROSS-VALIDATION
# USING 457 CORRELATION-SELECTED GENES
############################################################

# 1. set the working directory
setwd("~/Downloads/DepMap_Project") 

# 2. load the required packages
library(data.table)
library(randomForest)

# 3. load processed datasets
final_dataset <- fread("Dataset_A_gene_expression_ic50.csv")
X <- fread("filtered_gene_matrix.csv")

# 4. define target variable
y <- final_dataset$log_ic50

# 5. convert feature matrix to data frame
X <- as.data.frame(X)

# 6. check dimensions
cat("Number of samples =", nrow(X), "\n")
cat("Number of selected genes =", ncol(X), "\n")

# 7. create metric function
calculate_metrics <- function(actual, predicted) {
  
  mae <- mean(abs(predicted - actual))
  rmse <- sqrt(mean((predicted - actual)^2))
  r2 <- cor(predicted, actual)^2
  
  return(c(MAE = mae, RMSE = rmse, R2 = r2))
}

# 8. create 5 folds
set.seed(42)

k <- 5

folds <- sample(
  rep(1:k, length.out = nrow(X))
)

# 9. create empty data frames
cv_results <- data.frame(
  Fold = integer(),
  MAE = numeric(),
  RMSE = numeric(),
  R2 = numeric()
)

cv_predictions <- data.frame(
  Actual = numeric(),
  Predicted = numeric(),
  Fold = integer()
)

# 10. run 5-fold cross-validation
for (i in 1:k) {
  
  cat("Running fold", i, "\n")
  
  test_idx <- which(folds == i)
  train_idx <- setdiff(1:nrow(X), test_idx)
  
  X_train_cv <- X[train_idx, ]
  X_test_cv  <- X[test_idx, ]
  
  y_train_cv <- y[train_idx]
  y_test_cv  <- y[test_idx]
  
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
      Fold = i,
      MAE = as.numeric(metrics_cv["MAE"]),
      RMSE = as.numeric(metrics_cv["RMSE"]),
      R2 = as.numeric(metrics_cv["R2"])
    )
  )
  
  cv_predictions <- rbind(
    cv_predictions,
    data.frame(
      Actual = y_test_cv,
      Predicted = pred_cv,
      Fold = i
    )
  )
}

# 11. reset row names
rownames(cv_results) <- NULL
rownames(cv_predictions) <- NULL

# 12. print fold results
print(cv_results)

# 13. calculate summary
cv_summary <- data.frame(
  Metric = c("MAE", "RMSE", "R2"),
  Mean = c(
    mean(cv_results$MAE),
    mean(cv_results$RMSE),
    mean(cv_results$R2)
  ),
  SD = c(
    sd(cv_results$MAE),
    sd(cv_results$RMSE),
    sd(cv_results$R2)
  )
)

print(cv_summary)

# 14. calculate overall cross-validated performance
overall_cv_mae <- mean(
  abs(cv_predictions$Predicted - cv_predictions$Actual)
)

overall_cv_rmse <- sqrt(
  mean((cv_predictions$Predicted - cv_predictions$Actual)^2)
)

overall_cv_r2 <- cor(
  cv_predictions$Predicted,
  cv_predictions$Actual
)^2

cat("Overall CV MAE =", overall_cv_mae, "\n")
cat("Overall CV RMSE =", overall_cv_rmse, "\n")
cat("Overall CV R² =", overall_cv_r2, "\n")

# 15. save cross-validation results
write.csv(
  cv_results,
  "rf_5fold_cv_results.csv",
  row.names = FALSE
)

write.csv(
  cv_summary,
  "rf_5fold_cv_summary.csv",
  row.names = FALSE
)

write.csv(
  cv_predictions,
  "rf_5fold_cv_predictions.csv",
  row.names = FALSE
)

# 16. save actual vs predicted plot
png(
  "rf_5fold_cv_actual_vs_predicted.png",
  width = 2400,
  height = 2000,
  res = 300
)

plot(
  cv_predictions$Actual,
  cv_predictions$Predicted,
  pch = 19,
  xlab = "Actual log10(IC50)",
  ylab = "Predicted log10(IC50)",
  main = "Random Forest 5-Fold Cross-Validation\nActual vs Predicted"
)

abline(
  0,
  1,
  lwd = 2,
  lty = 2
)

dev.off()

# 17. save RMSE by fold plot
png(
  "rf_5fold_cv_RMSE_by_fold.png",
  width = 2200,
  height = 1600,
  res = 300
)

barplot(
  cv_results$RMSE,
  names.arg = cv_results$Fold,
  xlab = "Fold",
  ylab = "RMSE",
  main = "Random Forest 5-Fold Cross-Validation RMSE by Fold"
)

dev.off()

############################################################