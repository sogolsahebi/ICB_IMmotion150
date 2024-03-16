# Format_downloaded_data.R

# This script formats and cleans clinical and expression data.
# - Creates "CLIN.txt" dimension 589 x 32
# - Creates "EXPR.txt.gz" dimension 60308 x 325

library(tidyr)
library(dplyr)
library(readr)
library(data.table)

# Read the CSV file for clinical data
clin_path <- "~/BHK lab/ICB_IMmotion150/files/run_sample.csv"
clin <- read_csv(clin_path)

# Convert sample_attributes into separate columns and rename 'Patient.Identifier' to 'patient'
clin <- clin %>%
  separate_rows(sample_attributes, sep = ";") %>%
  separate(sample_attributes, into = c("Key", "Value"), sep = "=") %>%
  spread(Key, Value)
colnames(clin)[colnames(clin) == "clin$sample_alias"] <- "patient"

# Save formatted clinical data as CLIN.txt
clin_output_path <- "~/BHK lab/ICB_IMmotion150/files/CLIN.txt"
write.table(clin, clin_output_path, quote = FALSE, sep = "\t", row.names = FALSE)

# Read and format expression data
expr_path <- "~/BHK lab/ICB_IMmotion150/files/modified_transcript_tpms_all_samples.tsv"
expr <- as.data.frame(fread(expr_path))

# Extract gene name and remove duplicates
expr$gene_name <- sapply(strsplit(as.character(expr$target_id), "\\|"), `[`, 6)
expr <- expr[!duplicated(expr$gene_name), ]
rownames(expr) <- expr$gene_name
expr$target_id <- NULL
expr$gene_name <- NULL

# Sort the row names of 'expr'
expr <- expr[sort(rownames(expr)), ]

# Confirm all column values are numeric (output not shown)
all(sapply(expr, is.numeric)) == TRUE

# Save formatted expression data as EXPR.txt.gz
expr_output_path <- "~/BHK lab/ICB_IMmotion150/files/EXPR.txt.gz"
gz_conn <- gzfile(expr_output_path, "w")
write.table(expr, gz_conn, sep = "\t", row.names = TRUE, quote = FALSE)
close(gz_conn)
