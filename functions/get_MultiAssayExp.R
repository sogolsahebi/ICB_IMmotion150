library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringr)
#add.
library(dplyr)
library(readr)

#source_location <- "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main"
# source_location <- "~/Documents/GitHub/Pachyderm/PredictIO_ICB/ICB_Common"

get_MultiAssayExp <- function(study, input_dir, expr_with_counts_isoforms=FALSE){ 
  #source(paste(source_location, "code", "CNA_Functions.R", sep = "/"))
  #source(paste(source_location, "code", "SNV_Functions.R", sep = "/"))
  #source(paste(source_location, "code", "Create_SummarizedExp.R", sep = "/"))
  #data = read.csv( paste(source_location, "data", "DATASET_LOAD_INFO.csv", sep="/") , sep=";" , stringsAsFactors=FALSE )
  
  #add
  #study = "Wolf"
  #input_dir = "~/BHK lab/ICB_Wolf/files/"
  
  #add
  data <- read_delim("files/DATASET_LOAD_INFO.csv", delim = ";",quote = ";", escape_double = FALSE, trim_ws = TRUE)
  colnames(data) <- gsub("\"", "", names(data))
  data <- data %>% mutate(across(everything(), ~str_replace_all(.x, "\"", "")))
  data <- data[data$study == study, ]
  
  se_list <- Create_SummarizedExperiments( 
    input_dir=input_dir,
    study= data$study, 
    expr_bool= data$expr_bool, 
    snv_bool= data$snv_bool, 
    cna_bool= data$cna_bool, 
    cin_bool= data$cin_bool, 
    #TODO: this was added beacuse of the ICB_IMmotion150
    coverage= as.numeric(data$coverage), 
    indel_bool= data$indel_bool,
    expr_with_counts_isoforms=expr_with_counts_isoforms
  )
  
  #add 
  #se_list[["expr"]] <- S1
  
  cols <-list()
  for(assay_name in names(se_list)){
    cols[[assay_name]]<- data.frame(colData(se_list[[assay_name]]))
  }
  
  # Format and merge coldata
  allcols <- lapply(cols, function(col){
    return(rownames(col))
  })
  allcols <- unique(unlist(allcols))
  
  coldata <- NA
  for(col in cols){
    if(!is.data.frame(coldata)){
      coldata <- col
    }
    missing <- allcols[!allcols %in% rownames(coldata)]
    filtered <- col[rownames(col) %in% missing, ]
    if(nrow(filtered) > 0){
      coldata <- rbind(coldata, filtered)
    }
  }
  coldata <- coldata[order(rownames(coldata)), ]
   
  #saveRDS(IM, file = "output data/ICB_IMmotion150.rds", compress = "xz", compression_level = 9)
  #saveRDS(IM, file = "~/BHK lab/ICB/ICB_IMmotion150/output data/ICB_IMmotion150.rds")
   
  return(MultiAssayExperiment(experiments=se_list, colData=coldata))
}

