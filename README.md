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

## 📂 Step 2: Load Required Data Files

In this step, we load all the necessary data files for the eQTL analysis:

- **Gene expression matrices** for liver and testis (CSV files with genes as rows and samples as columns).
- **Genotype matrix** containing SNP data with SNP IDs as row names.
- **Covariate metadata** to adjust for confounding effects.
- **SNP and gene position files** which contain the genomic locations required by Matrix eQTL.

Each file should have samples or features as row names, depending on the data type.

```r
# Load gene expression data for liver
LIver_expr <- read.csv("expression_data_liver.csv", row.names = 1, check.names = FALSE)
LIver_expr <- data.matrix(LIver_expr)

# Load gene expression data for testis
Test_expr <- read.csv("expression_data_Testis.csv", row.names = 1, check.names = FALSE)
Test_expr <- data.matrix(Test_expr)

# Load genotype matrix (SNPs as rows, samples as columns)
geno <- read.csv("Genotype_matrix.csv", row.names = 1, check.names = FALSE)
geno <- data.matrix(geno)

# Simplify SNP names by removing allele suffixes (e.g., "_A")
geno_names <- rownames(geno)
rownames(geno) <- sub("_[ACGT]+$", "", geno_names)

# Load covariate metadata (samples as rows, covariates as columns)
cov <- read.csv("Metadata.csv", row.names = 1, check.names = FALSE)
cov <- data.matrix(cov)

# Extract second column as the covariate vector (e.g., batch or treatment)
cov_vector <- setNames(cov[,2], rownames(cov))

# Convert covariate vector to 1-row matrix for Matrix eQTL input
cov_matrix <- matrix(cov_vector, nrow = 1)

# Assign sample names as column names of the covariate matrix
colnames(cov_matrix) <- names(cov_vector)

# Load SNP position file (required by Matrix eQTL)
snpspos <- read.csv("SNP_Position_File.csv", header = TRUE, row.names = 1)

# Load gene position file (required by Matrix eQTL)
genepos <- read.csv("Gene_Position_file.csv", header = TRUE, row.names = 1)
```
---


---
## 🔄 Step 3: Convert Data to SlicedData Objects for Matrix eQTL

Matrix eQTL requires the input data to be in the form of **SlicedData** objects. This step converts your loaded expression, genotype, and covariate matrices into these objects.

```r
# Convert liver gene expression matrix to SlicedData object
gene_liver <- SlicedData$new()
gene_liver$CreateFromMatrix(LIver_expr)

# Convert testis gene expression matrix to SlicedData object
gene_testis <- SlicedData$new()
gene_testis$CreateFromMatrix(Test_expr)

# Convert genotype matrix to SlicedData object
snps <- SlicedData$new()
snps$CreateFromMatrix(geno)

# Convert covariate matrix to SlicedData object
cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(cov_matrix)
```
---

---

## ⚙️ Step 4: Run eQTL Analysis Using Matrix eQTL

In this step, we perform the actual eQTL mapping using the `Matrix_eQTL_main()` function. This function tests associations between SNP genotypes and gene expression levels, adjusting for covariates.

---

### What are cis and trans eQTLs?

- **Cis-eQTLs**: These are SNPs located near the gene they regulate, typically within 1 megabase (1 million base pairs) upstream or downstream. Cis-eQTLs often have direct regulatory effects on gene expression.
- **Trans-eQTLs**: These are SNPs that affect genes located far away on the genome, possibly on different chromosomes. Trans-eQTLs reflect distal regulatory relationships and are usually harder to detect due to weaker effect sizes and more multiple testing burden.

---

### Running cis-eQTL Analysis for Liver

```r
me_cis_Liver <- Matrix_eQTL_main(
  snps = snps,
  gene = gene_liver,
  cvrt = cvrt,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name = NULL,                    # no trans output here
  output_file_name.cis = "cis_eqtls_liver.txt",  # save cis-eQTL results
  pvOutputThreshold = 0,                      # skip trans tests
  pvOutputThreshold.cis = 1e-1,               # p-value threshold for cis (lenient for exploratory)
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 1e6,                              # Define cis-window as 1 megabase
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE
)
```
### Running trans-eQTL Analysis for Liver
```r
me_trans_Liver <- Matrix_eQTL_main(
  snps = snps,
  gene = gene_liver,
  cvrt = cvrt,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name = "trans_eqtls_liver.txt",  # save trans-eQTL results
  output_file_name.cis = NULL,                 # skip cis tests here
  pvOutputThreshold = 1e-6,                    # strict p-value threshold for trans (to control false positives)
  pvOutputThreshold.cis = 0,                   # disable cis tests
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 0,                                 # no cis window for trans
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE
)
```

### Running cis-eQTL Analysis for Testis

```r
me_cis_Testis <- Matrix_eQTL_main(
  snps = snps,
  gene = gene_testis,
  cvrt = cvrt,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name = NULL,                      # no trans output here
  output_file_name.cis = "cis_eqtls_testis.txt",  # save cis eQTL results
  pvOutputThreshold = 0,                        # skip trans
  pvOutputThreshold.cis = 1e-2,                 # enable cis p-value threshold
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 1e6,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE
)
```

### Running trans-eQTL Analysis for Testis

```r
me_trans_Testis <- Matrix_eQTL_main(
  snps = snps,
  gene = gene_testis,
  cvrt = cvrt,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name = "trans_eqtls_testis.txt",  # save trans eQTL results
  output_file_name.cis = NULL,                   # skip cis
  pvOutputThreshold = 1e-6,                      # enable trans p-value threshold
  pvOutputThreshold.cis = 0,                     # disable cis
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 0,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE
)
```

---

---
## 📊 Step 5: Extract and Visualize Significant eQTL Results

After running the eQTL analysis, we extract significant cis- and trans-eQTLs based on a p-value threshold and visualize their distributions using histograms and QQ plots. This helps to assess the overall significance and distribution of the associations.

---

### Extract and Filter cis-eQTLs

```r
# Extract cis-eQTL results
cis_eqtls_Liver <- me_cis_Liver$cis$eqtls
cis_eqtls_Liver <- as.data.table(cis_eqtls_Liver)
# Filter for significant cis-eQTLs (p-value < 0.01)
cis_eqtls_Liver <- cis_eqtls_Liver %>% filter(pvalue < 0.01)

# Extract trans-eQTL results
trans_eqtls_Liver <- me_trans_Liver$all$eqtls
trans_eqtls_Liver <- as.data.table(trans_eqtls_Liver)
# Filter for significant trans-eQTLs (p-value < 0.01)
trans_eqtls_Liver <- trans_eqtls_Liver %>% filter(pvalue < 0.01)


######## For Testis ##################
# Extract cis-eQTL results
cis_eqtls_Testis <- me_cis_Testis$cis$eqtls
cis_eqtls_Testis <- as.data.table(cis_eqtls_Testis)
# Filter for significant cis-eQTLs (p-value < 0.01)
cis_eqtls_Testis <- cis_eqtls_Testis %>% filter(pvalue < 0.01)

# For trans eQTLs
# Extract trans-eQTL results
trans_eqtls_Testis <- me_trans_Testis$all$eqtls
trans_eqtls_Testis <- as.data.table(trans_eqtls_Testis)
# Filter for significant trans-eQTLs (p-value < 0.01)
trans_eqtls_Testis <- trans_eqtls_Testis %>% filter(pvalue < 0.01)


# Plot histogram of cis/trans-eQTL p-values

par(mfrow = c(2, 2))

hist(cis_eqtls_Liver$pvalue, breaks=50, main="Histogram of cis-eQTL (Liver)", xlab="p-value")
hist(trans_eqtls_Liver$pvalue, breaks=100, main="Histogram of trans-eQTL (Liver)" , xlab="p-value")
hist(me_cis_Testis_$pvalue, breaks=50, main="Histogram of cis-eQTL (Testis)", xlab="p-value")
hist(trans_eqtls_Testis$pvalue, breaks=50, main="Histogram of trans-eQTL (Testis)", xlab="p-value")

```

- This QQ plot compares observed cis-eQTL p-values to expected p-values under no association. Points above the red line show more significant results than expected by chance, indicating true genetic effects. It helps check if the results are reliable and not due to random noise.
  
```r
qqplot(-log10(runif(nrow(cis_eqtls_Liver))), -log10(cis_eqtls_Liver$pvalue),
       main="QQ plot of cis-eQTL p-values",
       xlab="Expected -log10(p)", ylab="Observed -log10(p)")
abline(0,1, col="red")

```




