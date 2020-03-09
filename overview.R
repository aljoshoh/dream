library(org.Hs.eg.db)
library(dplyr)
library(ComplexHeatmap)

library(synapser)
library(synapserutils)
synLogin("aljoshoh","pw")

hs <- org.Hs.eg.db
ensemble <- 
  AnnotationDbi::select(hs, 
       keys = as.character(target_list$target),
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL") # TEC gene is doubled

rnaseq <- read.csv("dream_data/rnaseq.csv")
rna <- rnaseq
row.names(rna) = rnaseq$Gene
rna <- rna[,-c(1,2)] %>% t
dnaseq <- read.csv("dream_data/dnaseq.csv")
dnaseq$value <- 1
dnaseq <- spread(data = dnaseq[,c("lab_id","Hugo_Symbol","value")], key = Hugo_Symbol, value = value)
auc <- read.csv("dream_data/aucs.csv")
auc <- spread(data=auc, key=inhibitor, value=auc)
row.names(auc) = make.names(auc$lab_id)
auc <- auc[,-1]

clin_response <- read.csv("dream_data/response.csv")
clin_feat_cat <- read.csv("dream_data/clinical_categorical.csv")
clin_feat_legend <- read.csv("dream_data/clinical_categorical_legend.csv")
clin_feat_num <- read.csv("dream_data/clinical_numerical.csv")


tmp <- table(dnaseq$lab_id, dnaseq$Hugo_Symbol)
### Use string database to filter for informative genes by network propagation of drug targets
### Use morgan fingerprints of drugs that are used



View(list$score)


## PHONG
test <-  AnvSigGen(as.matrix(rna[,1:20]), as.matrix(auc))



test <- loadRData("features/drugs/drug_class.RData")







library(dplyr)
# dont forget make.names !
#### VALIDATION PHASE

rna_original <- read.csv("dream_data_leaderboard/rnaseq.csv",sep=",") #CHECK
rna <- read.csv("features_validation/rna/rnaseq_full.csv", sep=",")
rownames<- rna$Gene
rna <- rna[,-c(1,2)]
rna <- sapply(rna, as.numeric)
row.names(rna) = rownames
rna <- t(rna)

auc_original <- read.csv("dream_data_leaderboard/aucs.csv",sep=",")
auc <- read.csv("features_validation/auc/aucs_full.csv", sep=",")
auc <- auc[,-1] #!!not in original data <-------------------------------------------------------------
auc <- spread(data=auc, key=inhibitor, value=auc)
row.names(auc) = make.names(auc$lab_id)
auc <- auc[,-1]

mut_original <- read.csv("dream_data_leaderboard/dnaseq.csv",sep=",") #CHECK
mut <- read.csv("features_validation/mut/dnaseq_full.csv")
mut$value <- 1
mut <- mut[-c(576),]
mut <- spread(data = mut[,c("lab_id","var_name","value")], key = var_name, value = value)
mut[is.na(mut)] <- 0
row.names(mut) <- make.names(mut$lab_id)
mut <- mut[,-1]
mut <- mut[, apply(mut, 2, sum) >= 4 ]


clin_cat <- read.csv("dream_data_leaderboard/clinical_categorical.csv",sep=",") # Dream data
clin_cat <- clin_cat %>% mutate_all(factor)
row.names(clin_cat) <- make.names(clin_cat$lab_id)
clin_cat <- clin_cat[,-1]
clin_num <- read.csv("dream_data_leaderboard/clinical_numerical.csv")
row.names(clin_num) <- make.names(clin_num$lab_id)
clin_num <- clin_num[,-1]
clin <- read.csv("features_validation/clinical_data/full_clinical_data.csv") # Engineered training data
row.names(clin) = as.character(clin$b_id)
clin <- clin[,-1]
row.names.clin <- row.names(clin)
clin <- clin %>% mutate_at(colnames(clin_cat), factor)
row.names(clin) = make.names(row.names.clin)
clin$treatment_stratification <- factor(clin$treatment_stratification)
clin <- model.matrix(~., data=clin)

resp <- read.csv("features_validation/survival/response_full.csv",sep=",")
resp_original <- read.csv("dream_data_leaderboard/response.csv")
resp <- resp[,-1] # not in original
row.names(resp) = make.names(resp$lab_id)
resp <- resp[,-1]
resp$vitalStatus <- as.numeric(resp$vitalStatus)

save(clin, file = feature_path)
save(auc, file = response_path)

