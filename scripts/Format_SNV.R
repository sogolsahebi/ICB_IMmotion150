# Load the necessary library

# Format_SNV.R
# - Creates "SNV.csv" 

library(data.table)

# Set the path for the snv file and read the data using fread
snv_path <- "~/BHK lab/ICB_IMmotion150/files/COMBINED_snv.tsv.gz"
snv <- fread(snv_path, sep = "\t") # dimension 193692932 x 9

#Column names are already correct:("Sample", "Gene", "Chr", "Pos", "Ref", "Alt", "Effect", "MutType")]

# We are interested in keeping these Effects only according to line 201: https://github.com/BHKLAB-DataProcessing/ICB_Common/blob/main/code/Create_SummarizedExp.R
# Example output of unique(snv$Effect) is omitted for brevity

# Assuming 'snv' is already read as a data.table
# Make a copy of snv
snv_data <- copy(snv)
# Define categories for recoding 'Effect' values
splice_site_effects <- c("disruptive_inframe_deletion", "splice_acceptor_variant&intron_variant",
                         "splice_donor_variant&intron_variant", "splice_region_variant",
                         "splice_region_variant&intron_variant", "splice_region_variant&non_coding_transcript_exon_variant",
                         "splice_region_variant&synonymous_variant")

in_frame_del_effects <- c("frameshift_variant", "frameshift_variant&splice_region_variant",
                          "frameshift_variant&start_lost", "frameshift_variant&stop_gained")

in_frame_ins_effects <- c("conservative_inframe_insertion", "disruptive_inframe_insertion",
                          "disruptive_inframe_insertion&splice_region_variant")

missense_mutation_effects <- c("missense_variant", "missense_variant&splice_region_variant")

# Apply recodings using data.table's vectorized operations
snv_data[Effect %in% splice_site_effects, Effect := "Splice_Site"]
snv_data[Effect %in% in_frame_del_effects, Effect := "In_Frame_Del"]
snv_data[Effect %in% in_frame_ins_effects, Effect := "In_Frame_Ins"]
snv_data[Effect == "start_lost", Effect := "Start_Codon_Ins"]
snv_data[Effect %in% missense_mutation_effects, Effect := "Missense_Mutation"]
snv_data[Effect == "stop_lost", Effect := "Stop_Codon_Del"]
snv_data[Effect %in% "start_lost&conservative_inframe_deletion&splice_region_variant", Effect := "De_novo_Start_OutOfFrame"]
snv_data[Effect %in% "start_lost&splice_region_variant", Effect := "Start_Codon_SNP"]

# Write the modified SNV data to a file
path <- "~/BHK lab/ICB_IMmotion150/files/SNV.csv"
write.table(snv_data, path, quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)


# OPTIONAL RECOMMENDED: since the file is huge subseting the snv data to required effects values 
# We are interested in keeping these Effects only according to line 201: https://github.com/BHKLAB-DataProcessing/ICB_Common/blob/main/code/Create_SummarizedExp.R

desired_effects <- c("In_Frame_Del", "In_Frame_Ins", "Start_Codon_Ins", "Frame_Shift_Del",
                     "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                     "Splice_Site", "Stop_Codon_Del", "De_novo_Start_OutOfFrame", "Start_Codon_SNP")

# Filter the dataset to keep only rows with desired effects
snv_data <- snv_data[Effect %in% desired_effects] # now 7978686 x 9

path <- "~/BHK lab/ICB_IMmotion150/files/Subset_SNV.csv"
write.table(snv_data, path, quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)










