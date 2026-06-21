# Integrative Modeling of Doxorubicin Response in Lung Cancer Cell Lines

This repository contains the code, report, and poster for a senior design project focused on modeling doxorubicin response in lung cancer cell lines using gene expression and drug response data.

The project integrates **CCLE gene expression data** with **PRISM doxorubicin IC50 data** and applies machine learning methods to identify a candidate transcriptomic signature associated with doxorubicin sensitivity.

## Project Aim

The main objective of this project is to identify a gene expression-based predictive signature for doxorubicin response in lung cancer cell lines.

The project aims to:

* integrate CCLE gene expression data with PRISM doxorubicin IC50 data,
* filter the integrated dataset for lung cancer cell lines,
* transform IC50 values into log10(IC50),
* select genes associated with doxorubicin response,
* build Random Forest regression models,
* apply Recursive Feature Elimination (RFE),
* identify a final predictive gene signature,
* interpret the selected genes using functional enrichment analysis.

## Background

Doxorubicin is a widely used chemotherapeutic drug. However, cancer cells can show different levels of sensitivity or resistance to doxorubicin. One possible reason for this variation is molecular heterogeneity, including differences in gene expression profiles.

Gene expression data can provide information about the biological state of cancer cells. By combining transcriptomic data with drug response measurements, this project investigates whether specific gene expression patterns are associated with doxorubicin sensitivity or resistance.

In this project, doxorubicin response was represented using **IC50 values**. IC50 is the drug concentration required to inhibit cell viability by 50%. Higher IC50 values generally indicate lower sensitivity or higher resistance, while lower IC50 values indicate higher sensitivity.

## Data Sources

The datasets used in this project were obtained from the DepMap project:

* **CCLE gene expression dataset**
* **PRISM Repurposing Drug Screen dose-response dataset**
* **DepMap sample information dataset**

Raw datasets are not included in this repository due to file size limitations.

They can be downloaded from:

https://depmap.org

## Dataset Summary

After dataset integration and filtering, the final modeling dataset was constructed as follows:

| Processing Step                                               |         Output |
| ------------------------------------------------------------- | -------------: |
| PRISM doxorubicin entries with available IC50 values          | 471 cell lines |
| Common cell lines between CCLE expression and PRISM IC50 data | 460 cell lines |
| Lung cancer cell lines after filtering                        |  96 cell lines |
| Initial gene expression features                              |   19,220 genes |
| Genes selected by correlation filtering                       |      457 genes |
| Final selected predictive signature                           |       60 genes |

## Project Workflow

### 1. Dataset Integration

CCLE gene expression data and PRISM doxorubicin IC50 data were matched using the shared **DepMap ID**.

Main steps:

* load CCLE gene expression data,
* load PRISM drug response data,
* filter PRISM data for doxorubicin,
* remove samples with missing IC50 values,
* identify common DepMap IDs,
* merge gene expression data with IC50 values,
* filter for lung cancer cell lines,
* transform IC50 values into log10(IC50).

The final response variable used for modeling was **log10(IC50)**.

### 2. Correlation-Based Feature Selection

Since the original gene expression matrix contained thousands of genes, feature selection was performed before model training to reduce dimensionality.

Main steps:

* remove non-gene columns,
* remove zero-variance genes,
* calculate Pearson correlation between each gene and log10(IC50),
* select genes with absolute Pearson correlation values greater than or equal to 0.25.

This step reduced the feature set from **19,220 genes** to **457 genes**.

### 3. Baseline Random Forest Regression

A baseline Random Forest regression model was built using the 457 correlation-selected genes.

Model evaluation was performed using **5-fold cross-validation**.

Performance metrics included:

* MAE,
* RMSE,
* R².

In this project, R² was calculated as the squared Pearson correlation between actual and predicted log10(IC50) values.

### 4. Recursive Feature Elimination

Recursive Feature Elimination was applied to identify smaller and more informative gene sets.

The RFE procedure used Random Forest feature importance values based on **%IncMSE**. Genes with lower importance were removed iteratively, and multiple candidate gene sets were generated.

Each RFE-selected gene set was evaluated using **5-fold cross-validation**.

### 5. Final Model Selection

Two main criteria were considered:

* **Overall RMSE**, which measures prediction error,
* **Overall R²**, which measures explanatory performance.

The 8-gene model achieved the lowest Overall RMSE. However, the 60-gene model achieved the highest Overall R² and provided better biological interpretability.

Therefore, the **60-gene Random Forest RFE model** was selected as the final predictive signature.

## Main Results

| Model                    | Gene Set Size | Overall RMSE | Overall R² | Selection Reason                          |
| ------------------------ | ------------: | -----------: | ---------: | ----------------------------------------- |
| Baseline Random Forest   |     457 genes |       0.3027 |     0.2315 | Initial model after correlation filtering |
| Best RMSE RFE model      |       8 genes |       0.2773 |     0.3488 | Lowest prediction error                   |
| Final selected RFE model |      60 genes |       0.2790 |     0.3747 | Highest Overall R²                        |

The final 60-gene signature included genes such as:

* **GRN**
* **NUP98**
* **MAP2K2**
* **KLF4**
* **RBM39**
* **EEF2**
* **WIPI2**

These genes may be related to biological mechanisms such as survival signaling, apoptosis regulation, nuclear transport, MAPK signaling, RNA splicing, protein synthesis, and autophagy.

## Functional Enrichment Analysis

Functional enrichment analysis was performed using **g:Profiler** with **Homo sapiens** as the selected organism.

The final 60-gene signature was mainly associated with cellular component and protein complex terms, including:

* cytosol,
* organelle subcompartment,
* late endosome,
* trans-Golgi network transport vesicle membrane,
* granular vesicle,
* Grb2-Sos complex.

These findings suggest that doxorubicin response may involve multiple cellular processes rather than a single drug-response pathway.

## Repository Structure

```text
doxorubicin-response-modeling/
│
├── code/
│   ├── 01_dataset_integration_CCLE_PRISM.R
│   ├── 02_correlation_feature_selection.R
│   ├── 03_baseline_rf_5fold_cv.R
│   └── 04_rf_rfe_all_gene_sets_5fold_cv.R
│
├── poster/
│   └── senior_project_poster.pdf
│
├── report/
│   └── senior_project_report.pdf
│
└── README.md
```

## How to Run

The scripts should be run in the following order:

1. `01_dataset_integration_CCLE_PRISM.R`
2. `02_correlation_feature_selection.R`
3. `03_baseline_rf_5fold_cv.R`
4. `04_rf_rfe_all_gene_sets_5fold_cv.R`

Before running the scripts, download the required CCLE, PRISM, and sample information datasets from DepMap and place them in the working directory specified in the R scripts.

## Output Files

The scripts generate:

* integrated gene expression and IC50 dataset,
* filtered gene expression matrix,
* selected gene lists,
* baseline Random Forest cross-validation results,
* RFE gene set cross-validation summaries,
* final selected gene lists,
* model performance plots,
* actual versus predicted plots.

These outputs were used to prepare the final report and poster.

## Tools and Packages

This project was implemented in **R**.

Main tools and packages:

* R
* data.table
* randomForest
* Pearson correlation analysis
* Random Forest regression
* Recursive Feature Elimination
* 5-fold cross-validation
* g:Profiler

## Limitations

This study has several limitations:

* The final dataset included only 96 lung cancer cell lines.
* The model was developed using cell line data, which may not fully represent patient tumors.
* External validation was not performed.
* The selected 60-gene signature should be considered a candidate predictive signature.

## Future Work

Future work may include:

* validation using independent datasets,
* testing the model on patient-derived data,
* comparison with other machine learning algorithms,
* integration of additional omics data such as mutation or copy number profiles,
* experimental validation of selected genes.

## Authors

Cansu Akdemir\
Efehan Evren\
Esra Rençberler

Department of Genetics and Bioengineering\
Istanbul Bilgi University

## Project Information

Senior Design Project\
BIOE 492\
Department of Genetics and Bioengineering\
Istanbul Bilgi University

Supervised by Asst. Prof. Dr. Özlem Ulucan Açan

## Keywords

Doxorubicin, lung cancer, CCLE, PRISM, DepMap, IC50, gene expression, transcriptomics, Random Forest, Recursive Feature Elimination, pharmacogenomics, machine learning, drug response prediction.
