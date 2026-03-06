# Integrative Modeling of Doxorubicin Response in Cancer Cells

This project investigates how gene expression influences sensitivity to the chemotherapeutic drug doxorubicin.
The study integrates **gene expression–based drug response prediction** with **dynamic modeling of apoptosis signaling**.

The objective is to connect **molecular expression profiles** with **drug response metrics (IC50)** and explore how these responses translate into **apoptotic dynamics in cancer cells**.

---

# Background

Doxorubicin is a widely used chemotherapeutic agent that induces apoptosis by damaging DNA and activating the p53 signaling pathway. However, different cancer cell lines exhibit substantial variability in their sensitivity to this drug.

Gene expression profiles may help predict drug sensitivity across cell lines. At the same time, dynamic mathematical models can provide insight into how apoptosis signaling pathways respond to drug exposure.

This project combines **statistical modeling and systems biology approaches** to investigate doxorubicin response.

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

### 1. Dataset Construction

Gene expression and drug response datasets are integrated to construct a modeling dataset.

Steps:

* Load CCLE gene expression data
* Load PRISM drug response data
* Filter drug response data for **doxorubicin**
* Identify common cancer cell lines across datasets
* Merge gene expression and IC50 data
* Prepare the dataset for statistical modeling

Output:

A dataset containing gene expression features and corresponding **IC50 values** for doxorubicin.

---

### 2. Statistical Modeling of Drug Sensitivity

The constructed dataset is used to investigate the relationship between gene expression and drug sensitivity.

Methods explored include:

* linear regression modeling
* feature selection techniques
* identification of genes associated with drug response

The aim is to determine whether gene expression patterns can predict **IC50 values across cancer cell lines**.

---

### 3. Dynamic Modeling of Apoptosis

To connect drug sensitivity with biological mechanisms, a simplified dynamic model of the **p53–BAX–BCL2 apoptosis pathway** is explored.

Using ordinary differential equations (ODEs), the model investigates how drug-induced signaling may lead to apoptosis under different parameter conditions.

This step provides a **systems-level perspective on drug response**.

---

# Tools

* R
* data.table
* statistical modeling
* regression analysis
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
Project description and documentation

---

# Author

Cansu Akdemir

Genetics and Bioengineering / Istanbul Bilgi University

