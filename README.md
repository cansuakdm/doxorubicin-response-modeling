# Integrative Modeling of Doxorubicin Response in Cancer Cells

This project investigates how gene expression influences sensitivity to the chemotherapeutic drug doxorubicin.

The study focuses on gene expression–based drug response modeling to connect molecular data with drug sensitivity. The objective is to link gene expression profiles to drug response metrics (IC50) and evaluate whether these profiles can be used to predict response across cancer cell lines.

---

## Background

Doxorubicin is a widely used chemotherapeutic agent that induces apoptosis through DNA damage and activation of the p53 signaling pathway. However, cancer cell lines exhibit substantial variability in their sensitivity to this drug.

Gene expression data can be used to identify molecular features associated with drug response. Due to the high dimensionality of transcriptomic data, feature selection and regularized regression methods are required to build stable predictive models.

This project applies data-driven statistical modeling to investigate the relationship between gene expression and drug sensitivity.

---

## Data Sources

The datasets used in this project are obtained from the DepMap project:

- CCLE (Cancer Cell Line Encyclopedia) gene expression dataset  
- PRISM Repurposing Drug Screen dataset  

Due to file size limitations, raw datasets are not included in this repository.

They can be downloaded from:  
https://depmap.org

---

## Project Workflow

The project consists of three main stages.

---

### 1. Dataset Construction

Gene expression and drug response datasets are integrated to construct a unified modeling dataset.

Steps:

- Load CCLE gene expression data  
- Load PRISM drug response data  
- Filter for doxorubicin  
- Identify common cancer cell lines  
- Merge gene expression and IC50 data  
- Apply log transformation to IC50 values  

Output:

A dataset containing gene expression features and corresponding log(IC50) values.

---

### 2. Feature Selection and Statistical Analysis

The dataset is analyzed to identify genes associated with drug sensitivity.

Methods:

- Removal of zero-variance genes  
- Pearson correlation analysis between gene expression and log(IC50)  
- Threshold-based feature selection (|correlation| ≥ 0.3)  

This step reduces dimensionality and identifies candidate genes potentially associated with drug resistance and sensitivity.

Output:

A filtered gene expression matrix used for downstream modeling.

---

### 3. LASSO Regression Modeling

To predict drug response, a regularized regression model is applied.

Methods:

- Train-test split (80% training, 20% testing)  
- LASSO regression using glmnet  
- Cross-validation to determine optimal lambda  
- Prediction of log(IC50) values on the test set  

This step performs both prediction and automatic feature selection, selecting a subset of genes with non-zero coefficients.

---

## Tools

- R  
- data.table  
- glmnet  
- correlation analysis  
- feature selection  
- regression modeling  

---

## Repository Structure

code/  
R scripts used for dataset construction and analysis  

poster/  
Senior project poster  

report/  
Full senior project report  

README.md  
Project documentation  

---

## Author

Cansu Akdemir  

Genetics and Bioengineering  
Istanbul Bilgi University  
