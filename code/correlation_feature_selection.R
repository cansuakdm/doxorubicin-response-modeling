############################################################

# set the working directory
setwd("~/Downloads/DepMap_Project") 

# load the required package
library(data.table)

# load datasets
final_dataset <- fread("Dataset_A_gene_expression_ic50.csv")

# remove zero-variance genes
X <- final_dataset[, !c("depmap_id","ic50","log_ic50"), with=FALSE]
X <- X[, sapply(X, sd, na.rm=TRUE) > 0, with=FALSE]

# correlation
cor_results <- cor(X, final_dataset$log_ic50, method="pearson")

# convert correlation results to a vector
cor_vec <- cor_results[, 1]

# plot distribution of correlation values
png("gene_correlation_distribution.png", width = 800, height = 650)

hist(cor_vec,
     breaks = 50,
     main = "Distribution of Gene-log(IC50) Correlations",
     xlab = "Correlation",
     col = "lightblue",
     border = "black")

# add threshold lines (+0.1 and -0.1) and zero reference line
abline(v = 0.3, col = "darkgreen", lwd = 1.5, lty = 2)
abline(v = -0.3, col = "darkgreen", lwd = 1.5, lty = 2)
abline(v = 0, col = "red", lwd = 1.5)

dev.off()

# select genes with |correlation| >= 0.1 (threshold-based filtering)
selected_genes <- cor_vec[abs(cor_vec) >= 0.3]

# number of selected genes
length(selected_genes)

# create initial feature matrix using selected genes
X_filtered <- X[, names(selected_genes), with = FALSE]

# check dimensions
dim(X_filtered)

# save filtered matrix
fwrite(X_filtered, "filtered_gene_matrix.csv")

# save selected gene names
write.csv(names(selected_genes), "selected_genes.csv", row.names = FALSE)

############################################################
