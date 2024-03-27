# Format_downloaded_data.R

# This script formats and cleans clinical and expression data.
# - Creates "CLIN.txt" dimension 589 x 32
# - Creates 'expr_list.rds' including :
    # - EXPR_gene_tpm.tsv: dimension 46312 x 325
    # - EXPR_gene_counts.tsv: dimension 46312   325
    # - EXPR_tx_tpm.tsv: dimension 246624 x 325
    # - EXPR_tx_counts.tsv: dimension 246624 x 325

library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(tximport)

# Read the CSV file for clinical data
clin_path <- "~/BHK lab/ICB_IMmotion150/files/run_sample.csv"
clin <- read_csv(clin_path)

# Convert sample_attributes into separate columns and rename 'Patient.Identifier' to 'patient'
clin <- clin %>%
  separate_rows(sample_attributes, sep = ";") %>%
  separate(sample_attributes, into = c("Key", "Value"), sep = "=") %>%
  spread(Key, Value)
colnames(clin)[colnames(clin) == "sample_accession_id"] <- "patient"

# Save formatted clinical data as CLIN.txt
clin_output_path <- "~/BHK lab/ICB_IMmotion150/files/CLIN.txt"
write.table(clin, clin_output_path, quote = FALSE, sep = "\t", row.names = FALSE)

# Read and format expression data
# output: expr_list.rds
load("~/BHK lab/Annotation/Gencode.v19.annotation.RData")
work_dir <- "~/BHK lab/ICB_IMmotion150/files/kallisto_v0.46.1_GRCh38.40/"

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

unlink(file.path(work_dir, 'rnaseq', '__MACOSX'), recursive = TRUE)
process_kallisto_output(work_dir, tx2gene)

# SNV.txt.gz
snv <- read_excel(file.path(work_dir, '1-s2.0-S0092867417311224-mmc3.xlsx'), sheet='Table S3')
colnames(snv) <- snv[3, ]
snv <- snv[-c(1:3), ]
numcols <- c('Start', 'End', 'Tcov', 'Tac', 'Taf')
snv[, numcols] <- sapply(snv[, numcols], as.numeric)
gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep=";" , col.names=TRUE, row.names = FALSE )
close(gz)



