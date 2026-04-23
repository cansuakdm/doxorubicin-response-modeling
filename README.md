Integrative Modeling of Doxorubicin Response in Cancer Cell Lines

This project investigates how gene expression influences sensitivity to the chemotherapeutic drug doxorubicin in cancer cell lines.

The study focuses on building a predictive modeling pipeline that links gene expression profiles to drug response measurements, specifically IC50 and log(IC50) values. By integrating gene expression data with drug screening results, the project aims to identify molecular features associated with drug sensitivity and drug resistance, and to develop regression models that can predict response to doxorubicin.

The overall goal is to combine feature selection and machine learning–based regression modeling to better understand and predict variation in doxorubicin response across cancer cell lines.

Background

Doxorubicin is a widely used chemotherapeutic agent that exerts its effects primarily through DNA damage and apoptosis-related pathways. However, different cancer cell lines show substantial variability in their sensitivity to this drug.

Gene expression profiles provide a high-dimensional molecular snapshot of each cell line and can be used to identify biomarkers associated with treatment response. Because genome-scale expression datasets contain thousands of genes, feature selection is necessary to reduce dimensionality and retain the most informative predictors.

This project uses correlation-based filtering followed by regularized regression modeling to identify predictive genes and estimate drug response from gene expression data.

Data Sources

The datasets used in this project are obtained from the DepMap project:

CCLE (Cancer Cell Line Encyclopedia) gene expression dataset
PRISM Repurposing Drug Screen dataset containing doxorubicin response measurements

Due to file size limitations, raw datasets are not included in this repository.

They can be downloaded from:
https://depmap.org

Project Workflow

The project consists of three main stages.

1. Dataset Construction

Gene expression and drug response datasets are integrated to construct a unified modeling dataset.

Steps:

Load CCLE gene expression data
Load PRISM drug response data
Filter for doxorubicin
Identify common cancer cell lines
Merge gene expression data with IC50 values
Apply log transformation to IC50 values

Output:

A final dataset containing gene expression features together with IC50 and log(IC50) response variables.

2. Feature Selection and Statistical Analysis

Because the dataset contains thousands of genes, an initial filtering step is applied before modeling.

Methods:

Removal of zero-variance genes
Pearson correlation analysis between each gene and log(IC50)
Threshold-based feature selection using absolute correlation
Creation of a filtered gene matrix for downstream regression

In the current implementation, genes are selected using:

|correlation| ≥ 0.3

This step reduces dimensionality and identifies candidate genes potentially associated with doxorubicin sensitivity or resistance.

Output:

Filtered gene expression matrix
Selected gene list
Correlation distribution plot
3. Regression Modeling of Drug Response

The filtered gene set is used to build predictive models for log(IC50).

Methods:

Train/test split
Regularized regression using glmnet
Comparison of different model types:
LASSO (alpha = 1)
Elastic Net (alpha = 0.75, 0.5, 0.25)
Comparison of different numbers of top-ranked genes
Model evaluation using:
RMSE
R²

The best-performing configuration in the current analysis was:

Top 30 genes
Elastic Net with alpha = 0.25
Lambda = 0.04098384
Test R² = 0.493
RMSE = 0.276
24 selected genes in the final model

This stage identifies the most predictive subset of genes and provides a quantitative model of doxorubicin response.

Output:

Model comparison table
Best model coefficients
Best selected genes
Predicted vs. actual response plot
Tools
R
data.table
glmnet
Pearson correlation analysis
feature selection
regularized regression modeling
predictive model evaluation
Repository Structure

code/
R scripts used for dataset construction, feature selection, and regression modeling

poster/
Senior project poster

report/
Full senior project report

README.md
Project documentation

Current Findings

Preliminary results suggest that a subset of gene expression features can explain a meaningful proportion of the variation in doxorubicin response across cancer cell lines.

Among the tested models, an Elastic Net regression model outperformed standard LASSO settings in terms of test-set performance. This indicates that a partially regularized model that balances feature selection and coefficient shrinkage may better capture the gene-expression patterns associated with drug response.

These findings support the idea that transcriptomic data can be used to build predictive models for chemotherapy sensitivity.

Future Directions

Possible next steps include:

testing additional feature selection strategies
using repeated cross-validation for more robust performance estimates
comparing Elastic Net with other machine learning models
performing pathway-level interpretation of selected genes
validating the model on independent datasets
Author

Cansu Akdemir
Genetics and Bioengineering
Istanbul Bilgi University
