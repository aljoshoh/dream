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






#### VALIDATION PHASE
rna <- read.csv("features_validation/rna/rnaseq_full.csv", sep=",")
row.names(rna) <- rna$Gene
rna <- rna[,-c(1,2)]
rna <- t(rna)

auc <- read.csv("features_validation/auc/aucs_full.csv", sep=",")
auc <- auc[,-1]
auc <- spread(data=auc, key=inhibitor, value=auc)
row.names(auc) = make.names(auc$lab_id)
auc <- auc[,-1]

mut <- read.csv("features_validation/mut/dnaseq_full.csv")
mut$value <- 1
mut <- mut[-c(576),]
mut <- spread(data = mut[,c("lab_id","var_name","value")], key = var_name, value = value)
mut[is.na(mut)] <- 0
row.names(mut) <- make.names(mut$lab_id)
mut <- mut[,-1]
mut <- mut[, apply(mut, 2, sum) >= 4 ]



save(rna, file = feature_path)
save(auc, file = response_path)

