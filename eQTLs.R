


######### Eqtls analysis

library("MatrixEQTL")
library(qvalue)
library(GEOquery)
library(MatrixEQTL)
library(data.table)

library(data.table)

library(data.table)

# Path to folder with HTSeq count files
main_folder <- "/Eqtl/GSE93734_RAW"

# List only files that end with HTSeq_counts.txt.gz
count_files <- list.files(main_folder, pattern = "HTSeq_counts\\.txt\\.gz$", full.names = TRUE)

# Initialize expression list
expression_list <- list()
gene_ids <- NULL

# Loop through files
for (file in count_files) {
  # Extract sample name
  sample_name <- gsub("\\.HTSeq_counts\\.txt\\.gz$", "", basename(file))
  
  # Read data
  count_data <- fread(file, header = FALSE, col.names = c("GeneID", "Count"))
  
  # Remove special HTSeq rows (starting with "__")
  count_data <- count_data[!grepl("^__", GeneID)]
  
  # Store gene IDs once
  if (is.null(gene_ids)) {
    gene_ids <- count_data$GeneID
  }
  
  # Add expression counts
  expression_list[[sample_name]] <- count_data$Count
}

# Combine into matrix
expression_matrix <- as.data.table(expression_list)
expression_matrix <- cbind(GeneID = gene_ids, expression_matrix)

# Preview first few lines
print(head(expression_matrix))

# Save
fwrite(expression_matrix, "merged_expression_matrix.csv")




# Set file paths
expression_file <- "GSE113318_ProcessedMatrix.csv.gz"
covariates_file <- "GSE113318_Sample_Map.txt.gz"  # optional
genotypes <- fread("/Eqtl/GSE113318_for_R.raw")
# 1. Read expression data



# Transpose if needed: ensure genes are in rows, samples in columns
# If samples are in rows and genes in columns, do:
# expression_data <- transpose(expression_data)
expression_data =expression_matrix
# Remove the first column and filter by liver-related column names
# Ensure expression_data is a numeric matrix
expression_data <- as.matrix(expression_data)  # Ensure the data is numeric

# Remove the first column and filter by liver-related column names
expression_data_liver <- expression_data[, -1]  # Remove the first column
expression_data_liver <- expression_data_liver[, grepl("liver", colnames(expression_data_liver), ignore.case = TRUE)]  # Filter for liver columns
rownames(expression_data_liver)=rownames(expression_matrix)
# Check the result
head(expression_data_liver)
expression_data_liver=as.numeric(expression_data_liver)
# To view the complete matrix, you can use this:
dim(expression_data_liver)  # To check the dimensions of the filtered matrix

# Check the result (you can view the first few rows)
head(expression_data_liver)

# To view the complete matrix, you can use this:
dim(expression_data_liver)  # To check the dimensions of the filtered matrix
type(expression_data_liver)
# Check the result
head(expression_data_liver)ames and column headers
genes <- expression_data[[1]]
expression_matrix <- as.matrix(expression_data[,-1, with=FALSE])
rownames(expression_matrix) <- genes

# 2. Convert PLINK to dosage format using PLINK command line
# In terminal (not R), run:
# plink --bfile GSE113318_PLINK --recode A --out GSE113318_dosage

# 3. Read genotype dosage matrix

# Remove first 6 columns (FID, IID, etc.)
genotype_matrix <- as.matrix(genotypes[,-(1:6), with=FALSE])
rownames(genotype_matrix) <- genotypes$IID  # sample IDs as rownames
genotype_matrix <- t(genotype_matrix)  # SNPs in rows, samples in columns

# 4. Optional: Load covariates
# Covariates must be in format: rows = covariates, columns = samples
# If not using covariates, set modelCovariates = NULL
covariates <- NULL
# Example: if using covariates
# cov_data <- fread(covariates_file)
# Process and format as needed:
# covariates <- SlicedData$new()
# covariates$CreateFromMatrix(as.matrix(...))
clean_sample_ids <- sub(".*Sample_([0-9]+).*", "\\1", colnames(expression_data_liver))
colnames(expression_data_liver) <- clean_sample_ids


# 5. Prepare data for Matrix eQTL
expression_matrix=expression_data_liver

snps <- SlicedData$new()
# Find common sample IDs
common_samples <- intersect(colnames(expression_matrix), colnames(genotype_matrix))

# Subset both matrices to keep only those samples (columns)
expression_matrix_filtered <- expression_matrix[, common_samples, drop = FALSE]
genotype_matrix_filtered   <- genotype_matrix[, common_samples, drop = FALSE]











######################33 Cis trasn #############################################


library(MatrixEQTL)

# 3333333333333333################# adTA


###### Metadata 
GSE93734_metadata = getGEO("GSE93734",GSEMatrix = TRUE)
GSE93734_metadata$`GSE93734-GPL19176_series_matrix.txt.gz`
pheno_data_GSE93734 <- pData(GSE93734_metadata[[1]])

# Extract sample ID from the first column
pheno_data_GSE93734$SampleID <- sub("Sample ([0-9]+).*", "\\1", pheno_data_GSE93734$title)
# Extract the group (assuming it's already correct in the second column)
pheno_data_GSE93734$Group_Clean <- tolower(pheno_data_GSE93734$`group:ch1`)  # just to standardize
# View the cleaned metadata
clean_meta <- pheno_data_GSE93734[, c("SampleID", "group:ch1")]
# Get sample IDs from the gene expression matrix
expression_samples_colnames <- colnames(expression_matrix_filtered)
expression_samples_colnames=as.numeric(expression_samples_colnames)
clean_meta$SampleID=as.numeric(clean_meta$SampleID)
any(duplicated(filtered_meta$SampleID))
# Keep only metadata rows that are in the expression matrix
filtered_meta <- clean_meta[clean_meta$SampleID %in% expression_samples_colnames, ]

setdiff(clean_meta$SampleID, expression_samples_colnames)
setdiff(expression_samples_colnames, clean_meta$SampleID)

type(clean_meta$SampleID)
clean_meta$SampleID=unique(clean_meta$SampleID)

#######




library(rtracklayer)
gtf_sus=import("/Eqtl/Sus_scrofa.Sscrofa11.1.113.gtf")
snp_map <- read.table("/Eqtl/GSE113318_SNP_Map.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(snp_map)
bim <- fread("/Eqtl/GSE113318_binary.bim")

# Format: snp, chr, pos
snp_pos_df <- data.frame(
  snp = bim$V2,
  chr = bim$V1,
  pos = bim$V4
)

genes_gtf <- gtf_sus[gtf_sus$type == "gene"]

gene_pos_df <- data.frame(
  geneid = mcols(genes_gtf)$gene_id,
  chr = as.character(seqnames(genes_gtf)),
  start = start(genes_gtf),
  end = end(genes_gtf)
)

snp_pos_df <- data.frame(
  snpid = snp_map$Name,                # SNP ID (e.g., ALGA0000009)
  chr = snp_map$Chromosome,            # Chromosome
  pos = snp_map$Position            # Position (numeric)
                # Just a placeholder, can be removed if unused
)
# gene_pos_df$end <- gene_pos_df$start + 1000 
gene_pos_df <- gene_pos_df[, c("geneid", "chr", "start", "end")]

snp_pos_df <- snp_map[, c("Name", "Chromosome", "Position")]
colnames(snp_pos_df) <- c("snpid", "chr", "pos")

snp_pos_df$chr <- as.character(snp_pos_df$chr)

genotype_snp_ids <- sub("_[ACGT]+$", "", rownames(genotype_matrix_filtered))
# Add those to a new column so we can compare
snp_row_map <- data.frame(
  full_id = rownames(genotype_matrix_filtered),
  snp_id = genotype_snp_ids
)

# Keep only SNPs that are present in snp_pos_df$snpid
snp_row_map_filtered <- snp_row_map[snp_row_map$snp_id %in% snp_pos_df$snpid, ]

# Subset genotype matrix to those SNPs
genotype_matrix_filtered <- genotype_matrix_filtered[snp_row_map_filtered$full_id, ]


# Optional: subset snp_pos_df to match
snp_pos_df <- snp_pos_df[snp_pos_df$snpid %in% snp_row_map_filtered$snp_id, ]



library(biomaRt)
ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")

gene_pos_df <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = rownames(expression_matrix),
  mart = ensembl
)

colnames(gene_pos_df) <- c("geneid", "chr", "start", "end")


library(MatrixEQTL)

# SNPs
SNP_file <- SlicedData$new()
SNP_file$CreateFromMatrix(as.matrix(genotype_matrix_filtered))
rownames(SNP_file) <- sub("_[ACGT]+$", "", rownames(SNP_file))


# Expression
expression_matrix=expression_data_liver
expression_matrix_clean <- as.matrix(trimws(expression_matrix))
# Convert the matrix to numeric values
expression_matrix_clean <- matrix(as.numeric(expression_matrix_clean), nrow = nrow(expression_matrix), ncol = ncol(expression_matrix))
# Check if the conversion worked
head(expression_matrix_clean)
rownames(expression_matrix_clean) <- rownames(expression_matrix)


gene_ids_in_expression <- rownames(expression_matrix_clean)
gene_ids_in_pos <- gene_pos_df$geneid
# Find any mismatches
mismatched_genes <- setdiff(gene_ids_in_expression, gene_ids_in_pos)
if(length(mismatched_genes) > 0) {
  cat("Mismatched gene IDs:", mismatched_genes, "\n")
}
common_genes <- intersect(rownames(expression_matrix_clean), gene_pos_df$geneid)
expression_matrix_filtered <- expression_matrix_clean[common_genes, , drop = FALSE]
# Filter gene position data to only include common genes
gene_pos_filtered <- gene_pos_df[gene_pos_df$geneid %in% common_genes, ]


# Create the SlicedData object from the filtered expression matrix
gene_file <- SlicedData$new()

gene_file$CreateFromMatrix(expression_matrix_filtered)



# Set parameters based on the paper
cisDist <- 1e6  # 1 Mb cutoff for cis-acting eQTLs
pvThreshold_cis <- 1e-3  # P-value threshold for cis-acting eQTLs
pvThreshold_trans <- 1e-6  # P-value threshold for trans-acting eQTLs



# Proceed with eQTL analysis
me <- Matrix_eQTL_main(
  snps = SNP_file,
  gene = gene_file,
  cvrt = cvrt,
  output_file_name = NULL,
  pvOutputThreshold = 1e-6,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  snpspos = snp_pos_df,
  genepos = gene_pos_filtered,
  cisDist = 1e6,
  pvOutputThreshold.cis = 1e-3,
  verbose = TRUE
)

cis_eqtls <- me$cis$eqtls
cis_eqtls_significant <- cis_eqtls[cis_eqtls$FDR < 1e-6, ]

trans_eqtls <- me$trans$eqtls
trans_eqtls_significant <- trans_eqtls[trans_eqtls$FDR < 1e-6, ]





ggplot(cis_eqtls_significant, aes(x = snps, y = -log10(pvalue))) +
  geom_point() +
  theme_minimal() +
  xlab("SNP") +
  ylab("-log10(p-value)") +
  ggtitle("Cis-eQTLs Significant Associations")






#############


######################



library(vcfR)
vcf <- read.vcfR("/Eqtl/GSE113318_AGES_genotype.vcf")
vcf@meta
# Extract genotypes from the VCF file
genotypes <- vcf@gt  # Genotypes are stored in the "gt" slot

# Calculate genotyping rate (proportion of non-missing genotypes per SNP)
genotyping_rate <- rowMeans(is.na(genotypes))  # Proportion of missing values for each SNP

# Filter SNPs with genotyping rate < 5% (remove SNPs with more than 5% missing data)
vcf_filtered_genotyping_rate <- vcf[genotyping_rate < 0.05, ]


genotypes_numeric <- apply(genotypes, 1, function(x) {
  # Convert "0/0" -> 0, "0/1" -> 1, "1/1" -> 2, etc.
  sapply(strsplit(x, "/"), function(y) ifelse(y[1] == y[2], as.numeric(y[1]), 1))
})

# Function to perform Hardy-Weinberg equilibrium test

# Function to perform Hardy-Weinberg equilibrium test
# Function to perform Hardy-Weinberg equilibrium test
# Example function for Hardy-Weinberg equilibrium test
HWE_test <- function(snp) {
  # Skip SNPs with only one unique genotype (monomorphic SNP)
  if(length(unique(snp)) == 1) {
    return(NA)  # Return NA for monomorphic SNPs (cannot be tested)
  }
  
  # Run Hardy-Weinberg equilibrium test
  p_value <- tryCatch({
    # HWE test for the current SNP (use your preferred method or package)
    # Here, we assume genotypes are coded as 0, 1, 2 for alleles
    chisq.test(table(factor(snp, levels = c(0, 1, 2))))$p.value
  }, error = function(e) {
    return(NA)  # Return NA in case of an error
  })
  
  return(p_value)  # Return the p-value from HWE test
}



# Apply HWE test to all SNPs
hwe_pvalues <- apply(genotypes, 1, HWE_test)

# Filter SNPs based on HWE p-value threshold (e.g., 1e-6)
vcf_filtered_hwe <- vcf_filtered_genotyping_rate[hwe_pvalues > 1e-6, ]

# Step 3: MAF Filter
maf <- apply(genotypes, 1, function(x) min(table(x) / length(x)))
vcf_filtered_maf <- vcf_filtered_hwe[maf > 0.05, ]

# Step 4: Save the Filtered VCF File
write.vcf(vcf_filtered_maf, "filtered_genotype_data.vcf")




