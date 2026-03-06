# Integrative Modeling of Doxorubicin Response in Cancer Cells

This project investigates the relationship between gene expression and drug sensitivity in cancer cell lines. The goal is to integrate gene expression–based IC50 prediction with dynamic modeling of apoptosis.

## Background

Doxorubicin is a widely used chemotherapeutic drug that induces apoptosis in cancer cells. However, drug response varies significantly across cell lines. Gene expression profiles may help predict drug sensitivity, while dynamic models can provide insights into apoptosis signaling dynamics.

## Data Sources

The datasets used in this project are obtained from the DepMap project:

- CCLE gene expression dataset
- PRISM drug response dataset

Due to file size limitations, the raw datasets are not included in this repository.

They can be downloaded from:  
https://depmap.org

## Workflow

1. Load gene expression data (CCLE)
2. Load drug response data (PRISM)
3. Filter for doxorubicin
4. Identify common cell lines
5. Merge gene expression and IC50 data
6. Prepare dataset for statistical modeling

## Tools

- R
- data.table
- statistical modeling
- basic systems biology modeling

## Project Structure

code/ – R scripts for dataset creation and analysis  
poster/ – senior project poster  
README.md – project description

## Author

Cansu Akdemir  
Genetics and Bioengineering  
Istanbul Bilgi University
