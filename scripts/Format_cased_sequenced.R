# File: Format_cased_sequnced.R
# Goal: Save cased_sequenced.csv (dimensions:  589x  4 ).

# Parse command line arguments for input, output, and annotation directory paths.
#args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]

# Load the clinical merged data from the specified file path.
#clin <- read.table(file.path(work_dir, 'CLIN.txt'), sep="\t", header=TRUE)
clin_path <- "~/BHK lab/ICB_IMmotion150//files/CLIN.txt"
clin <- read.csv(clin_path, sep="\t", header=TRUE)

# sort the 'patient' column from clin and store it as Pateint.
patient <- sort(clin$patient)

# Initialize a data frame for 'case' with the unique patients and default values
case <- as.data.frame(cbind(patient, rep(0, length(patient)), rep(0, length(patient)), rep(0, length(patient))))
colnames(case) <- c("patient", "snv", "cna", "expr")
rownames(case) <- patient

# Convert the case values to numeric.
case$snv <- as.numeric(as.character(case$snv))
case$cna <- as.numeric(as.character(case$cna))
case$expr <- as.numeric(as.character(case$expr))

# Load the RNA data
# Read the RNA-Seq data from the gct file.
#expr <- read.csv(file.path(input_dir, "EXPR.txt.gz"), stringsAsFactors=FALSE , sep="\t" )
expr_path <- "~/BHK lab/ICB_IMmotion150//files/EXPR.txt.gz"
expr <- read.csv(expr_path, stringsAsFactors=FALSE , sep="\t",check.names = FALSE )

# Sort the row names of 'expr'
expr <- expr[sort(rownames(expr)),]

# Check the overlap of patient IDs between the 'expr' and 'clin' data
sum(colnames(expr) %in% clin$patient)  #325 overlap

# Update the 'expr' column in 'case' based on the presence of patient IDs in the 'expr' data
for(i in 1:nrow(case)) {
  if(rownames(case)[i] %in% colnames(expr)) {
    case$expr[i] = 1
  }
}

#TODO: add the SNV

# Save the updated 'case' data frame to a CSV file.
#write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
case_path <- "~/BHK lab/ICB_IMmotion150//files/cased_sequenced.csv"
write.table( case , case_path , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
