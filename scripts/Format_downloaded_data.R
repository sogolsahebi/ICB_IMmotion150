# Format_downloaded_data.R

# This script formats and cleans clinical and expression data.
# - Creates "CLIN.txt" dimension 589 x 32
# - Creates 'expr_list.rds' including :
    # - EXPR_gene_tpm.csv: dimension 46312 x 325
    # - EXPR_gene_counts.csv: dimension 46312   325
    # - EXPR_tx_tpm.csv: dimension 246624 x 325
    # - EXPR_tx_counts.csv: dimension 246624 x 325

library(tidyr)
library(dplyr)
library(data.table)
library(tximport)
library(data.table)

# Read the CSV file for clinical data
clin <- read.csv("files/run_sample.csv")

# Convert sample_attributes into separate columns and rename 'Patient.Identifier' to 'patient'
clin <- clin %>%
  separate_rows(sample_attributes, sep = ";") %>%
  separate(sample_attributes, into = c("Key", "Value"), sep = "=") %>%
  spread(Key, Value)
colnames(clin)[colnames(clin) == "sample_accession_id"] <- "patient"

# Save formatted clinical data as CLIN.txt
write.table(clin, "files/CLIN.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Read and format expression data
# output: expr_list.rds
load("~/BHK lab/Annotation/Gencode.v40.annotation.RData")
work_dir <- "files/kallisto_v0.46.1_GRCh38.40/"

dir.create(file.path(work_dir, 'rnaseq'))

# Creates EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv 
process_kallisto_output <- function(work_dir, tx2gene){
  
  samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names=FALSE, recursive = FALSE)
  files <- file.path(work_dir, 'rnaseq', samples, "abundance.h5")
  names(files) <- samples
  
  expr_tx <- tximport(files, type = "kallisto", txOut = TRUE, ignoreAfterBar = TRUE)
  expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  
  expr_list <- list()
  expr_list[['expr_gene_tpm']] <- log2(expr_gene$abundance + 0.001)
  expr_list[['expr_gene_counts']] <- log2(expr_gene$counts + 1)
  expr_list[['expr_isoform_tpm']] <- log2(expr_tx$abundance + 0.001)
  expr_list[['expr_isoform_counts']] <- log2(expr_tx$counts + 1)
  
  saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))
  
  unlink(file.path(work_dir, 'rnaseq'), recursive = TRUE)
}
process_kallisto_output(work_dir, tx2gene)

# SNV.txt.gz
snv <- fread("files/COMBINED_snv.tsv.gz", sep = "\t") # dimension 193692932 x 9

