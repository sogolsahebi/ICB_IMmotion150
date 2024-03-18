# Format_downloaded_data.R

# This script formats and cleans clinical and expression data.
# - Creates "CLIN.txt" dimension 589 x 32
# - Creates "EXPR.txt.gz" dimension 61544 x 326

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
colnames(clin)[colnames(clin) == "sample_alias"] <- "patient"

# Save formatted clinical data as CLIN.txt
clin_output_path <- "~/BHK lab/ICB_IMmotion150/files/CLIN.txt"
write.table(clin, clin_output_path, quote = FALSE, sep = "\t", row.names = FALSE)

# Read and format expression data
expr_path <- "~/BHK lab/ICB_IMmotion150/files/modified_transcript_tpms_all_samples.tsv"
expr <- as.data.frame(fread(expr_path))

# Extract gene_ids
expr$gene_id <- sapply(strsplit(as.character(expr$target_id), "\\|"), `[`, 2)

# Aggregating expression data to gene-level, focusing on numeric columns
# - Group by 'gene_id' to handle each gene as a unique entity
# - Summarise across all numeric columns to aggregate transcript-level data into gene-level
# - 'where(is.numeric)': Ensures only numeric columns (e.g., expression values) are considered for summing
# - 'na.rm = TRUE': Ignores NA values, preventing them from affecting the summation
# - Ungroup to remove the grouping structure, returning a regular data frame
expr_aggregated <- expr %>%
  group_by(gene_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# Convert tibble to a traditional data frame
expr_aggregated <- as.data.frame(expr_aggregated)

rownames(expr_aggregated) <- expr_aggregated$gene_id
expr_aggregated$gene_id <- NULL
expr_aggregated$target_id <- NULL

# Sort the row names of 'expr'
expr_aggregated <- expr_aggregated[sort(rownames(expr_aggregated)), ] #dim 61544 x 326

# Confirm all column values are numeric (output not shown)
all(sapply(expr_aggregated, is.numeric)) == TRUE

# Save formatted expression data as EXPR.txt.gz
expr_output_path <- "~/BHK lab/ICB_IMmotion150/files/EXPR.txt.gz"
gz_conn <- gzfile(expr_output_path, "w")
write.table(expr_aggregated, gz_conn, sep = "\t", row.names = TRUE, quote = FALSE)
close(gz_conn)
