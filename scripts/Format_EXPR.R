# Format_EXPR.R.

# RNA-seq Data Processing.
# - save "Format_EXPR.R "(dimensions 61544   325).

library(data.table)


# Data Reading
# Read EXPR.txt.gz file.
#expr <- read.csv(file.path(input_dir, "EXPR.txt.gz"), stringsAsFactors=FALSE , sep="\t" )
expr_path <- "~/BHK lab/ICB_IMmotion150/files/EXPR.txt.gz"
expr <- read.csv(expr_path, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

# Data Filtering
# Read 'cased_sequenced.csv' file
#case = read.csv( file.path(input_dir, "cased_sequenced.csv") , sep=";" )
file_path <- "~/BHK lab/Ravi_Testing/files/cased_sequenced.csv"
case <- read.csv(file_path, sep = ";")

# Filter the 'expr' dataset to include only patients with expr value of 1 in the 'case' dataset
expr <- expr[ , case[case$expr %in% 1,]$patient]

# Data Transformation
# Convert TPM data to log2-TPM for consistency with other data formats
expr <- log2(expr + 0.001) # range is now -9.965784 16.801317

# Save expr as csv file, EXPR.csv.
#write.table( tpm , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
file_path <- "~/BHK lab/ICB_IMmotion150/files/EXPR.csv"
write.table( expr , file_path , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
