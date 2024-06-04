library(data.table)

Get_TMB_raw = function( file , case ){
  #TODO: improvment dont need to read snveach time not efficient
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE)
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] , ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s)
  }
  tmb
}

Get_nsTMB_raw = function( file , case ){
  #snv = read.csv(  file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] & 
               snv$Effect %in% c("In_Frame_Del" , "In_Frame_Ins" , "Start_Codon_Ins" , "Frame_Shift_Del" ,
                                 "Frame_Shift_Ins" , "Missense_Mutation" , "Nonsense_Mutation" , "Nonstop_Mutation" ,
                                 "Splice_Site" , "Stop_Codon_Del" , "De_novo_Start_OutOfFrame" , "Start_Codon_SNP") ,
    ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s)
  }
  tmb
}

Get_indel_TMB_raw = function( file , case , indel_bool ){
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ snv$Sample %in% snv_patient[i] & snv$MutType %in% "INDEL" , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s)
    }
    
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_indel_nsTMB_raw = function( file , case  , indel_bool ){
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ 	snv$Sample %in% snv_patient[i] & 
                  snv$MutType %in% "INDEL" & 
                  snv$Effect %in% c("In_Frame_Del","In_Frame_Ins","Start_Codon_Ins","Frame_Shift_Del",
                                    "Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                    "Splice_Site","Stop_Codon_Del") , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s) 
    }
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_TMB_perMb = function( file , case , coverage ){
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] , ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
  }
  tmb	
}

Get_nsTMB_perMb = function( file , case , coverage ){
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] & 
               snv$Effect %in% c("In_Frame_Del","In_Frame_Ins","Start_Codon_Ins","Frame_Shift_Del",
                                 "Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                 "Splice_Site","Stop_Codon_Del") ,
    ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
  }
  tmb
}

Get_indel_TMB_perMb = function( file , case , coverage , indel_bool ){
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ snv$Sample %in% snv_patient[i] & snv$MutType %in% "INDEL" , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
    }
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_indel_nsTMB_perMb = function( file , case , coverage , indel_bool ){
  #snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  snv = file
  cases = read.csv( case , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ 	snv$Sample %in% snv_patient[i] & 
                  snv$MutType %in% "INDEL" & 
                  snv$Effect %in% c("In_Frame_Del","In_Frame_Ins","Start_Codon_Ins","Frame_Shift_Del",
                                    "Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                    "Splice_Site","Stop_Codon_Del") , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
    }
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_SNV_feature = function( case, file, coverage , indel_bool , mutsig_bool ){
  # case = case_file
  # file = snv_file
  # tmb_raw <- Get_TMB_raw( file=file , case=case )
  # nstmb_raw <- Get_nsTMB_raw( file=file , case=case )
  # indel_tmb_raw <- Get_indel_TMB_raw( file=file , case=case , indel_bool=indel_bool )
  # indel_nstmb_raw <- Get_indel_nsTMB_raw( file=file , case=case , indel_bool=indel_bool )
  # tmb_permb <- Get_TMB_perMb( file=file , case=case , coverage=coverage )
  # nstmb_permb <- Get_nsTMB_perMb( file=file , case=case , coverage=coverage)
  # indel_tmb_permb <- Get_indel_TMB_perMb( file=file , case=case , coverage=coverage , indel_bool=indel_bool )
  # indel_nstmb_permb <- Get_indel_nsTMB_perMb( file=file , case=case , coverage=coverage , indel_bool=indel_bool )
  x = cbind( 
    Get_TMB_raw( file=file , case=case ) ,
    Get_nsTMB_raw( file=file , case=case ) ,
    Get_indel_TMB_raw( file=file , case=case , indel_bool=indel_bool ) ,
    Get_indel_nsTMB_raw( file=file , case=case , indel_bool=indel_bool ) ,
    Get_TMB_perMb( file=file , case=case , coverage=coverage ) ,
    Get_nsTMB_perMb( file=file , case=case , coverage=coverage) ,
    Get_indel_TMB_perMb( file=file , case=case , coverage=coverage , indel_bool=indel_bool ) ,
    Get_indel_nsTMB_perMb( file=file , case=case , coverage=coverage , indel_bool=indel_bool )
  )
  colnames(x) = c( "TMB_raw" , "nsTMB_raw" , "indel_TMB_raw" , "indel_nsTMB_raw" , 
                   "TMB_perMb" , "nsTMB_perMb" , "indel_TMB_perMb" , "indel_nsTMB_perMb" )
  return(x)
}

