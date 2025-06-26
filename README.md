## 🧬 eQTL Analysis in R Using Matrix eQTL (Liver & Testis)

This tutorial provides a step-by-step guide to performing expression quantitative trait loci (eQTL) analysis using the [Matrix eQTL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339296/) R package. We use RNA-Seq expression data and SNP genotype data from two tissues: **liver** and **testis**.

---

## 🔧 Step 1: Install and Load Required Packages

We begin by installing and loading the necessary R libraries. If you haven’t installed the required packages, uncomment the install commands.

```r
# Install packages if not already installed
# BiocManager::install("org.Ss.eg.db")

library(qvalue)
library(GEOquery)
library(MatrixEQTL)
library(data.table)
library(qqman)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(ggvenn)
library(clusterProfiler)
library(org.Ss.eg.db)
library(biomaRt)
```

---

## 📂 Step 2: Load Gene Expression Data

We load the gene expression matrices for liver and testis. Each file should be in CSV format:

* Rows = Genes
* Columns = Samples
* First column = row names (gene IDs)

```r
# Read liver expression matrix
LIver_expr <- read.csv("expression_data_liver.csv", row.names = 1, check.names = FALSE)
LIver_expr <- data.matrix(LIver_expr)

# Read testis expression matrix
Test_expr <- read.csv("expression_data_Testis.csv", row.names = 1, check.names = FALSE)
Test_expr <- data.matrix(Test_expr)
```

---

## 🧬 Step 3: Load Genotype Data

The genotype matrix should:

* Have SNPs as rows
* Samples as columns
* SNP IDs as row names (e.g., rs1234\_A)

We also clean the SNP names to keep only the base SNP ID:

```r
# Load genotype matrix
geno <- read.csv("Genotype_matrix.csv", row.names = 1, check.names = FALSE)
geno <- data.matrix(geno)

# Simplify SNP names by removing allele suffix
geno_names <- rownames(geno)
rownames(geno) <- sub("_[ACGT]+$", "", geno_names)
```

---

## 🧾 Step 4: Load Covariate Metadata

Covariates (e.g., **batch**, **RIN**, **sex**) are important to adjust for potential confounding effects in eQTL analysis. Including these in the model helps to account for variation that is not due to genotype.

Your `Metadata.csv` file should:

* ✅ Have **sample IDs as row names**
* ✅ Contain one or more **columns of covariates**

In this step, we extract the **second column** as a covariate of interest (e.g., batch or treatment group) and convert it into a **1-row matrix** required by the Matrix eQTL format.

```r
# Load metadata containing covariates (e.g., batch, RIN)
cov <- read.csv("Metadata.csv", row.names = 1, check.names = FALSE)
cov <- data.matrix(cov)

# Extract second column as covariate vector (e.g., batch information)
cov_vector <- setNames(cov[,2], rownames(cov))

# Convert vector to 1-row matrix for Matrix eQTL
cov_matrix <- matrix(cov_vector, nrow = 1)

# Set column names to match sample IDs
colnames(cov_matrix) <- names(cov_vector)
```

✅ **Important**: Make sure that the **sample names in your covariate matrix** match exactly with the column names of your **expression** and **genotype** matrices.

➡️ You can also extract **multiple covariates** by adjusting the column selection:

```r
cov_matrix <- t(cov[, c(2, 3)])  # for two covariates
```
