# Clinical Data Processing
# Goal: save CLIN.csv (dimensions: 589 x 36).

# Libraries and Source Files
# Access necessary functions from ICB_common codes.
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

# Read library tibble for add_column function
library(tibble)

# Parse command line arguments for input, output, and annotation directory paths.
# args <- commandArgs(trailingOnly = TRUE)
# input_dir <- args[1]
# output_dir <- args[2]
# annot_dir <- args[3]

# Define the path to the source data
file_path <- "~/BHK lab/ICB_IMmotion150/files/CLIN.txt"
# Read CLIN.txt file
clin <- read.csv(file_path, stringsAsFactors=FALSE, sep="\t", dec=',')

# Set RNA-related columns based on 'RNA_All' values
clin$rna <- ifelse(is.na(clin$AnnoRNASampleID), NA, "rnaseq")
clin$rna_info <- ifelse(is.na(clin$AnnoRNASampleID), NA, "tpm")

# Set DNA-related columns based on 'DNA_All' values
clin$dna <- ifelse(is.na(clin$AnnoNormalWESID), NA, "wes")

# TODO: what should be in dna_info?
clin$dna_info <- ifelse(is.na(clin$AnnoNormalWESID), NA, "--")

# Select the required columns for further analysis
selected_cols <- c("patient", "Stage", "ARM", "BestResponse", "PFS", "rna", "rna_info", "dna", "dna_info")

# Combine selected columns with additional columns
clin <- cbind(clin[, selected_cols], "Renal Cell Carcinoma", NA, NA, 
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Reorder columns.
colnames(clin) <- c("patient", "stage", "drug_type", "recist", "pfs", "rna", "rna_info", 
                    "dna", "dna_info", "primary", "sex", "age", "histo", "response.other.info", "response", 
                    "t.pfs", "t.os", "os")

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

# Use the format_clin_data function for further formatting.
clin <- format_clin_data(clin, "patient", selected_cols, clin)


# Read 'curation_tissue.csv' file
# annotation_tissue <- read.csv(file.path(annot_dir, 'curation_tissue.csv'))
path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv"
annotation_tissue <- read.csv(path)

# TODO: double check curation_tissue.csv row: Mmotion150, Renal Cell Carcinoma, Kidney

# Annotate 'clin' using the 'annotation_tissue' data; Set tissueid column (survival_unit and survival_type columns are added in this step).
clin <- annotate_tissue(clin=clin, study='IMmotion150', annotation_tissue=annotation_tissue, check_histo=FALSE) 


# TODO: Double check with Sisira drug_types.
# "Atezo+Bev" "Sunitinib" "Atezo" --> Atezolizumab + Bevacizumab   Sunitinib Atezolizumab

# Set treatmentid after tissueid column, based on curation_drug.csv file
annotation_drug <- read.csv("https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_drug.csv")
clin <- add_column(clin, treatmentid=annotate_drug('IMmotion150', clin$drug_type, annotation_drug), .after='tissueid')


# TODO: double check Farnoosh: Sunitinib ---> targeted, Atezolizumab + Bevacizumab ---> "IO+combo", Atezolizumab ---> "PD-1/PD-L1"
# Sunitinib is a targeted therapy.
# Atezolizumab + Bevacizumab is a combination of immunotherapy and targeted therapy.
# Atezolizumab is an immunotherapy.
# Six categories of:  PD-1/PD-L1, CLA4 , IO+combo, IO+chemo, Chemo+targeted, targeted

# Set drug_type based on treatmentid
clin$drug_type[clin$treatmentid == "Sunitinib"] <- 'targeted'
clin$drug_type[clin$treatmentid == "Atezolizumab + Bevacizumab"] <- 'IO+combo'
clin$drug_type[clin$treatmentid == "Atezolizumab"] <- 'PD-1/PD-L1'

# Replace empty string values with NA
clin[clin == ""] <- NA

# Save the processed data as CLIN.csv file
file_path <- "~/BHK lab/Ravi_Testing/files/CLIN.csv"
write.table(clin, file=file_path, quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)
