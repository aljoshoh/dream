# Contains functions for importing the test data
library(Hmisc)
library(dplyr)
library(tidyr)
#....


get_features <- function(
  ### Loads features
  ######################################################
  path = NULL, # RData file of feature matrix, e.g. "features/features.RData"
  matrixfy = T
){
  path <- as.character(path)
  message(paste0('Loading ', path))
  if(matrixfy){
    features <- as.matrix(loadRData(path))
  } else {
    features <- loadRData(path)
  }
  
  return(features)
}######################################################


import_aucs <- function(
  path = "dream_data_leaderboard/aucs.csv"
){
  auc <- read.csv(path, sep=",")
  auc <- spread(data=auc, key=inhibitor, value=auc)
  row.names(auc) = make.names(auc$lab_id)
  auc <- auc[,-1]
  auc[] <- auc %>% mutate_all(function(x) impute(x))
  selected_features <- c("17-AAG (Tanespimycin)", "A-674563", "Afatinib (BIBW-2992)", 
                         "Alisertib (MLN8237)", "AT7519", "Axitinib (AG-013736)", "AZD1480", 
                         "Barasertib (AZD1152-HQPA)", "BEZ235", "BMS-345541", "Bortezomib (Velcade)", 
                         "Bosutinib (SKI-606)", "Cabozantinib", "Canertinib (CI-1033)", 
                         "Cediranib (AZD2171)", "CHIR-99021", "CI-1040 (PD184352)", "Crenolanib", 
                         "Crizotinib (PF-2341066)", "CYT387", "Dasatinib", "Doramapimod (BIRB 796)", 
                         "Dovitinib (CHIR-258)", "Elesclomol", "Erlotinib", "Flavopiridol", 
                         "Foretinib (XL880)", "GDC-0879", "GDC-0941", "Gefitinib", "GSK-1838705A", 
                         "GSK-1904529A", "GSK690693", "GW-2580", "Ibrutinib (PCI-32765)", 
                         "Idelalisib", "Imatinib", "INK-128", "JAK Inhibitor I", "JNJ-28312141", 
                         "JNJ-38877605", "JNJ-7706621", "KI20227", "KU-55933", "KW-2449", 
                         "Lapatinib", "Linifanib (ABT-869)", "LY-333531", "Masitinib (AB-1010)", 
                         "MGCD-265", "Midostaurin", "MK-2206", "MLN120B", "MLN8054", "Motesanib (AMG-706)", 
                         "Neratinib (HKI-272)", "NF-kB Activation Inhibitor", "Nilotinib", 
                         "NVP-ADW742", "NVP-TAE684", "Pazopanib (GW786034)", "PD173955", 
                         "Pelitinib (EKB-569)", "PHA-665752", "PHT-427", "PI-103", "Ponatinib (AP24534)", 
                         "PP242", "PRT062607", "Quizartinib (AC220)", "RAF265 (CHIR-265)", 
                         "Rapamycin", "Regorafenib (BAY 73-4506)", "Roscovitine (CYC-202)", 
                         "Ruxolitinib (INCB018424)", "S31-201", "Saracatinib (AZD0530)", 
                         "SB-431542", "Selumetinib (AZD6244)", "SGX-523", "SNS-032 (BMS-387032)", 
                         "Sorafenib", "STO609", "SU11274", "Sunitinib", "TG100-115", "Tivozanib (AV-951)", 
                         "Tofacitinib (CP-690550)", "Tozasertib (VX-680)", "Trametinib (GSK1120212)", 
                         "Vandetanib (ZD6474)", "Vargetef", "Vatalanib (PTK787)", "Vismodegib (GDC-0449)", 
                         "VX-745", "XAV-939", "YM-155")
  auc <- auc[,selected_features]
  return(auc)
}


import_rnaseq <- function(
  path = "dream_data_leaderboard/rnaseq.csv"
){
  rna <- read.csv(path,sep=",")
  rownames<- rna$Gene
  rna <- rna[,-c(1,2)]
  rna <- sapply(rna, as.numeric)
  row.names(rna) = rownames
  rna <- t(rna)
  return(rna)
}

import_dnaseq <- function(
  path = "dream_data_leaderboard/dnaseq.csv"
){
  selected_features <- c("ADAMTS7_p.H1024R", "AGTPBP1_p.X97_splice", "ASXL1_p.G646Wfs*12", 
                         "DDX60L_p.P1169Q", "DNMT3A_p.R882C", "DNMT3A_p.R882H", "ENAH_p.Q150K", 
                         "FLT3_p.D835H", "FLT3_p.D835Y", "IDH1_p.R132C", "IDH1_p.R132H", 
                         "IDH2_p.R140Q", "IDH2_p.R172K", "IMPG1_p.T226K", "JAK2_p.V617F", 
                         "KRAS_p.G13D", "LRRCC1_p.A6V", "MYO5B_p.D670V", "NACAD_p.S654L", 
                         "NPM1_p.W288Cfs*12", "NRAS_p.G12D", "NRAS_p.G13D", "NRAS_p.Q61H", 
                         "NRAS_p.Q61K", "SF3B1_p.K666N", "SLC39A12_p.G426V", "SP4_p.E7K", 
                         "SRSF2_p.P95H", "SRSF2_p.P95L", "TMPRSS6_p.R446W", "TPRX1_p.S200P", 
                         "TRIO_p.E2226G", "TRIO_p.S2221G", "TRPM3_p.V75F", "TRPM3_p.X1180_splice", 
                         "U2AF1_p.S34F", "WNK1_p.Q796H", "ZNF195_p.H152L", "ZNF285_p.L24V", 
                         "ZNF711_p.C318F")
  mut <- read.csv(path,sep=",")
  mut$value <- 1
  ###mut <- mut[-c(576),] 
  mut <- pivot_wider(data = mut[,c("lab_id","var_name","value")], names_from = var_name, values_from = value, names_repair = 'unique')
  #mut <- spread(data = mut[,c("lab_id","var_name","value")], key = var_name, value = value)
  mut <- as.data.frame(mut)
  rownames <- mut$lab_id
  mut[is.null(mut)] <- 0
  mut[is.na(mut)] <- 0
  row.names(mut) <- make.names(mut$lab_id)
  #mut <- mut[,-1]
  mut <- as.data.frame(lapply(mut, function(x) unique(unlist(x))))#
  not_in_data <- selected_features[!selected_features %in% colnames(mut)]
  adding <- as.data.frame(matrix(0, ncol = length(not_in_data), nrow = nrow(mut)))
  colnames(adding) = not_in_data
  mut <- cbind(mut, adding)
  row.names(mut) = rownames
  mut <- mut[,selected_features] %>% as.matrix()
  return(mut)
}

import_clin <- function( #TODO: check if dfs need to be imputed
  path_num = "dream_data_leaderboard/clinical_numerical.csv",
  path_cat = "dream_data_leaderboard/clinical_categorical.csv"
){
  clin_cat <- read.csv(path_cat,sep=",")
  clin_cat <- clin_cat %>% mutate_all(factor)
  row.names(clin_cat) <- make.names(clin_cat$lab_id)
  clin_cat <- clin_cat[,-1]
  clin_num <- read.csv(path_num,sep=",")
  row.names(clin_num) <- make.names(clin_num$lab_id)
  clin_num <- clin_num[,-1]
  clin <- cbind(clin_num, clin_cat)
  clin <- model.matrix(~., data=clin)
  return(clin)  
}

import_response <- function(
  path = "dream_data_leaderboard/response.csv"
){
  resp <- read.csv(path,sep=",")
  row.names(resp) = make.names(resp$lab_id)
  resp <- resp[,-1]
  resp$vitalStatus <- as.character(resp$vitalStatus)
  resp$vitalStatus[resp$vitalStatus == "Alive"] <- 1
  resp$vitalStatus[resp$vitalStatus == "Dead"] <- 0
  resp$vitalStatus <- as.numeric(resp$vitalStatus)
  return(resp)
}


convert_path <- function(
  path = NULL, # a single path where the dnn model is stored, e.g. /storage/groups/cbm01/workspace/dream_aml/submission/sc1/models/h2o_models/
  path.in.object = NULL # path from object, e.g. MODEL$param[[i]][[1]][[1]]
){
  res <- gsub("/storage/groups/cbm01/workspace/dream_aml/metadata/h2odnn/", path, path.in.object)
  return(res)
}


