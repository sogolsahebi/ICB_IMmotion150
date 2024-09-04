# Clinical Data Processing
# Goal: save CLIN.csv (dimensions: 589 x 49).

# Libraries and Source Files.
# Access necessary functions from ICB_common codes.
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

# Read library tibble for add_column function
library(tibble)

# Read CLIN.txt file
clin_orginal <- read.csv("files/CLIN.txt", stringsAsFactors=FALSE, sep="\t", dec=',')

# Open case to set rna and dna relate dcolumns
case <- read.csv( "files/cased_sequenced.csv" , sep=";" , stringsAsFactors=FALSE )
rna_patients <- case$patient[case$expr == 1]
dna_patients <- case$patient[case$snv == 1]

# Update clinical data for RNA sequencing patients
clin_orginal$rna[clin_orginal$patient %in% rna_patients] <- "rnaseq"
clin_orginal$rna_info[clin_orginal$patient %in% rna_patients] <- "tpm"

# Update clinical data for DNA sequencing patients
clin_orginal$dna[clin_orginal$patient %in% dna_patients] <- "wes"
clin_orginal$dna_info[clin_orginal$patient %in% dna_patients] <- "snv"

# Select the required columns for further analysis
selected_cols <- c("patient", "Stage", "ARM", "BestResponse", "PFS", "rna", "rna_info", "dna", "dna_info")

# Combine selected columns with additional columns
clin <- cbind(clin_orginal[, selected_cols], "Kidney", NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Reorder columns
colnames(clin) <- c("patient", "stage", "drug_type", "recist", "t.pfs", "rna", "rna_info", 
                    "dna", "dna_info", "primary", "sex", "age", "histo", "response.other.info", "response", 
                    "pfs", "t.os", "os")

# Modify stage values, remove "STAGE II" to "II" example.
clin$stage <- gsub("STAGE ", "", clin$stage)

# Calculate the response using Get_Response function.
clin$response = Get_Response(data = clin)

# Reorder columns.
clin <- clin[, c(
  "patient", "sex", "age", "primary", "histo", "stage", 
  "response.other.info", "recist", "response", "drug_type", "dna", "dna_info", "rna", "rna_info", "t.pfs", 
  "pfs", "t.os", "os"
)]

# Set the pfs NA
clin_orginal$pfs <- NA
clin_orginal$t.pfs <- NA

#TODO: set survival_unit as NA

# Use the format_clin_data function for further formatting.
clin <- format_clin_data(clin_orginal, "patient", selected_cols, clin)

# Read 'curation_tissue.csv' file
# annotation_tissue <- read.csv(file.path(annot_dir, 'curation_tissue.csv'))
path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv"
annotation_tissue <- read.csv(path)

# Annotate 'clin' using the 'annotation_tissue' data; Set tissueid column (survival_unit and survival_type columns are added in this step).
clin <- annotate_tissue(clin=clin, study='IMmotion150', annotation_tissue=annotation_tissue, check_histo=FALSE) 

# Set treatmentid after tissueid column, based on curation_drug.csv file
annotation_drug <- read.csv("https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_drug.csv")
clin <- add_column(clin, treatmentid=annotate_drug('IMmotion150', clin$drug_type, annotation_drug), .after='tissueid')

# Sunitinib ---> targeted, Atezolizumab + Bevacizumab ---> "IO+targeted", Atezolizumab ---> "PD-1/PD-L1"
# Atezolizumab + Bevacizumab is a combination of immunotherapy and targeted therapy.
# Six categories of:  PD-1/PD-L1, CLA4 , IO+combo, IO+chemo, Chemo+targeted, targeted ,IO+targeted.

# Set drug_type based on treatmentid
clin$drug_type[clin$treatmentid == "Sunitinib"] <- 'targeted'
clin$drug_type[clin$treatmentid == "Atezolizumab + Bevacizumab"] <- 'IO+targeted'
clin$drug_type[clin$treatmentid == "Atezolizumab"] <- 'PD-1/PD-L1'

# Replace empty string values with NA
clin[clin == ""] <- NA

# Save the processed data as CLIN.csv file
write.table(clin, file="files/CLIN.csv", quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)

