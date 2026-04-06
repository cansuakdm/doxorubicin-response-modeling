# Integrative Modeling of Doxorubicin Response in Cancer Cells

This project investigates how gene expression influences sensitivity to the chemotherapeutic drug doxorubicin.

The study integrates **gene expression–based drug response modeling** with **dynamic modeling of apoptosis signaling** to connect molecular data with biological mechanisms.

The objective is to link **gene expression profiles** to **drug response metrics (IC50)** and explore how these responses may translate into **apoptotic dynamics in cancer cells**.

---

# Background

Doxorubicin is a widely used chemotherapeutic agent that induces apoptosis through DNA damage and activation of the p53 signaling pathway. However, cancer cell lines exhibit substantial variability in their sensitivity to this drug.

Gene expression data can be used to identify molecular features associated with drug response. In addition, dynamic mathematical models can provide insight into how apoptosis signaling pathways respond to drug exposure.

This project combines **data-driven statistical analysis** with **systems biology modeling** to investigate doxorubicin response.

---

# Data Sources

The datasets used in this project are obtained from the DepMap project:

* **CCLE (Cancer Cell Line Encyclopedia)** gene expression dataset
* **PRISM Repurposing Drug Screen** dataset containing drug response metrics

Due to file size limitations, raw datasets are not included in this repository.

They can be downloaded from:
https://depmap.org

---

# Project Workflow

The project consists of three main stages.

---

## 1. Dataset Construction

Gene expression and drug response datasets are integrated to construct a unified modeling dataset.

Steps:

* Load CCLE gene expression data
* Load PRISM drug response data
* Filter for **doxorubicin**
* Identify common cancer cell lines
* Merge gene expression and IC50 data
* Apply log transformation to IC50 values

Output:

A dataset containing gene expression features and corresponding **log(IC50)** values.

---

## 2. Feature Selection and Statistical Analysis

The dataset is analyzed to identify genes associated with drug sensitivity.

Methods:

* Removal of zero-variance genes
* Pearson correlation analysis between gene expression and log(IC50)
* Threshold-based feature selection (|correlation| ≥ 0.1)

This step reduces dimensionality and identifies candidate genes potentially associated with **drug resistance and sensitivity**.

Output:

A filtered gene expression matrix used for downstream modeling.

---

## 3. Dynamic Modeling of Apoptosis

To provide a biological interpretation of drug response, a simplified dynamic model of the **p53–BAX–BCL2 apoptosis pathway** is explored.

Using ordinary differential equations (ODEs), the model investigates how drug-induced signaling may lead to apoptosis under different conditions.

This step provides a **mechanistic perspective** linking drug response to cellular behavior.

---

# Tools

* R
* data.table
* correlation analysis
* feature selection
* regression modeling
* basic systems biology modeling

---

# Repository Structure

code/  
R scripts used for dataset construction and analysis  

poster/  
Senior project poster  

report/  
Full senior project report  

README.md  
Project documentation  

---

# Author

Cansu Akdemir  

Genetics and Bioengineering  
Istanbul Bilgi University
