# Format_EXPR.R
# Goal: EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv 
#   - EXPR_gene_tpm.csv: dimension 46312 x 325
#   - EXPR_gene_counts.csv: dimension 46312 x 325
#   - EXPR_tx_tpm.csv: dimension 246624 x 325
#   - EXPR_tx_counts.csv: dimension 246624 x 325

# load libraries
library(data.table)
library(stringr)

# Define the input and output to dir
input_dir <- "~/BHK lab/ICB_IMmotion150/files/kallisto_v0.46.1_GRCh38.40/"
output_dir <- "~/BHK lab/ICB_IMmotion150/files/"

# extract expr_list including 4 expr tsv files
expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))

# Save all 4 expr as csv file
for(assay_name in names(expr_list)){
  write.table( 
    expr_list[[assay_name]], 
    file= file.path(output_dir, paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')), 
    quote=FALSE, 
    sep=";", 
    col.names=TRUE, 
    row.names=TRUE 
  )
}
