# eQTL Analysis in R Using Matrix eQTL with Liver and Testis Data

This tutorial walks through the steps for performing eQTL analysis using the [Matrix eQTL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339296/) R package. We use gene expression and genotype data from liver and testis tissue.

---

## 🔧 Step 1: Install and Load Required Packages

```r
# Uncomment to install if not already installed
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



# Liver expression matrix
LIver_expr <- read.csv("expression_data_liver.csv", row.names = 1, check.names = FALSE)
LIver_expr <- data.matrix(LIver_expr)

# Testis expression matrix
Test_expr <- read.csv("expression_data_Testis.csv", row.names = 1, check.names = FALSE)
Test_expr <- data.matrix(Test_expr)


geno <- read.csv("Genotype_matrix.csv", row.names = 1, check.names = FALSE)
geno <- data.matrix(geno)

# Fix SNP names if needed
rownames(geno) <- sub("_[ACGT]+$", "", rownames(geno))


cov <- read.csv("Metadata.csv", row.names = 1, check.names = FALSE)
cov <- data.matrix(cov)

# Convert to matrix format for Matrix eQTL
cov_vector <- setNames(cov[,2], cov[,1])
cov_matrix <- matrix(cov_vector, nrow = 1)

# Assign sample names
colnames(cov_matrix) <- names(cov_vector)


gene_liver <- SlicedData$new()
gene_liver$CreateFromMatrix(LIver_expr)

gene_testis <- SlicedData$new()
gene_testis$CreateFromMatrix(Test_expr)

snps <- SlicedData$new()
snps$CreateFromMatrix(geno)

cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(cov_matrix)



snpspos <- read.csv("SNP_Position_File.csv", header = TRUE, row.names = 1)
genepos <- read.csv("Gene_Position_file.csv", header = TRUE, row.names = 1)
