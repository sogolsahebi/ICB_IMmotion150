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
library(data.table)

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
load("~/BHK lab/Annotation/Gencode.v40.annotation.RData")
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
process_kallisto_output(work_dir, tx2gene)


# SNV.txt.gz

snv_path <- "~/BHK lab/ICB_IMmotion150/files/COMBINED_snv.tsv.gz"
snv <- fread(snv_path, sep = "\t") # dimension 193692932 x 9

# # Define the exact mapping for some effects
# effect_mapping <- c(
#   "disruptive_inframe_deletion" = "In_Frame_Del",
#   "conservative_inframe_insertion" = "In_Frame_Ins",
#   "start_lost" = "Start_Codon_Ins",
#   "frameshift_variant" = "Frame_Shift_Del", # Assume deletions unless specified as insertion
#   "frameshift_variant&stop_gained" = "Frame_Shift_Ins", # Specific case where it's clearly an insertion
#   "missense_variant" = "Missense_Mutation",
#   "stop_gained" = "Nonsense_Mutation",
#   "stop_lost" = "Nonstop_Mutation",
#   "splice_region_variant" = "Splice_Site",
#   "splice_donor_variant&splice_region_variant&intron_variant" = "Splice_Site",
#   "stop_lost&splice_region_variant" = "Stop_Codon_Del",
#   "start_lost&conservative_inframe_deletion&splice_region_variant" = "De_novo_Start_OutOfFrame",
#   "start_lost&splice_region_variant" = "Start_Codon_SNP"
# )
# 
# # Function to map the effects with a default behavior for unspecified frameshift variants
# map_effects <- function(effect) {
#   if (!is.na(effect_mapping[effect])) {
#     return(effect_mapping[effect])
#   } else if (grepl("frameshift_variant", effect)) {
#     # Additional check for frameshift variants to decide if it is a deletion or not specified
#     if (grepl("insertion", effect, ignore.case = TRUE)) {
#       return("Frame_Shift_Ins")
#     } else {
#       return("Frame_Shift_Del")  # Default to Frame_Shift_Del if not specified as insertion
#     }
#   }
#   return(effect)  # Return the effect unchanged if it does not match any specific mapping
# }
# 
# # Apply the mapping function to each effect in the dataframe
# snv$Effect <- sapply(snv$Effect, map_effects)
# 
# 
# # c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" )
# 
# # Print updated effects
# unique(snv$Effect)
# 
# subset_snv <- snv[ snv$Effect %in% c("In_Frame_Del" , "In_Frame_Ins" , "Start_Codon_Ins" , "Frame_Shift_Del" ,
#                                           "Frame_Shift_Ins" , "Missense_Mutation" , "Nonsense_Mutation" , "Nonstop_Mutation" ,
# #                                           "Splice_Site" , "Stop_Codon_Del" , "De_novo_Start_OutOfFrame" , "Start_Codon_SNP") ,]
# gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
# write.table( snv , file=gz , quote=FALSE , sep=";" , col.names=TRUE, row.names = FALSE )
# close(gz)




