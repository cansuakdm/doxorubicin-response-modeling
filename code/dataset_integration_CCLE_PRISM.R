############################################################

# 1. set the working directory
setwd("~/Downloads/DepMap_Project")   # change accordingly
getwd()

# 2. load the required package
library(data.table)

# 3. load datasets
expr <- fread("CCLE_expression.csv")
prism <- fread("prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv")

# 4. filter doxorubicin
doxo <- prism[name == "doxorubicin"]    # drug names are lowercase

# remove missing IC50 values
doxo <- doxo[!is.na(ic50)]

# check number of cell lines
nrow(doxo)

# 5. prepare expression dataset (rename first column to depmap_id)
setnames(expr, "V1", "depmap_id")

# 6. identify common cell lines
common_ids <- intersect(expr$depmap_id, doxo$depmap_id)
length(common_ids)

# 7. subset the expression dataset to keep only cell lines that have matching IC50 data
expr_subset <- expr[
  depmap_id %in% common_ids
]

# 8. subset IC50 data
doxo_subset <- doxo[
  depmap_id %in% common_ids,
  .(depmap_id, ic50)
]

# 9. merge expression and IC50
final_dataset <- merge(
  expr_subset,
  doxo_subset,
  by = "depmap_id"
)

# 9.5 load sample metadata (contains cancer type info)
sample_info <- fread("sample_info.csv")

# Check column names to identify cancer type column
colnames(sample_info)

# 9.6 filter lung cancer cell lines using primary_disease
lung_ids <- sample_info[
  grepl("lung", primary_disease, ignore.case = TRUE),
  DepMap_ID
]

# 9.7 keep only lung cancer samples in the final dataset
final_dataset <- final_dataset[
  depmap_id %in% lung_ids
]

# 9.8 check dataset dimensions after filtering
dim(final_dataset)
head(final_dataset)

# 10. convert IC50 to numeric (IC50 may be stored as character → convert to numeric)
final_dataset$ic50 <- as.numeric(final_dataset$ic50)

# 11. Log-transform IC50 (Log transformation improves regression modeling)
final_dataset$log_ic50 <- log10(final_dataset$ic50)

# 12. save Dataset A
fwrite(final_dataset, "Dataset_A_gene_expression_ic50.csv")

############################################################
