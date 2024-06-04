# File: Format_cased_sequnced.R
# Goal: Save cased_sequenced.csv (dimensions:  589x  4 ).

# Load the clinical merged data 
clin <- read.csv("files/CLIN.txt", sep="\t", header=TRUE)

# sort the 'patient' column from clin and store it as Patient.
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
expr_list <- readRDS(file.path("files/kallisto_v0.46.1_GRCh38.40/", 'expr_list.rds'))
expr <- expr_list[['expr_gene_tpm']]
rna_patients <- sort(colnames(expr))

# Load SNV data 
snv <- fread("files/COMBINED_snv.tsv.gz", sep = "\t") # dimension 193692932 x 9
dna_patient <- unique(snv$Sample)

# Check the overlap of patient IDs between the 'expr' and 'clin' data
sum(colnames(expr) %in% clin$patient)  #325 overlap 

# Check the overlap of patient IDs between the 'dna' and 'clin' data
sum(dna_patient %in% clin$patient) # 326 overal

# Update the 'expr' column in 'case' based on the presence of patient IDs in the 'expr' data
for(i in 1:nrow(case)) {
  if(rownames(case)[i] %in% rna_patients) {
    case$expr[i] = 1
  }
  if(rownames(case)[i] %in% dna_patient ) {
    case$snv[i] = 1
  }
}

write.table( case , "files/cased_sequenced.csv" , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
