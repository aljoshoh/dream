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
auc <- auc[,apply(auc,2, function(x) length(which(is.na(x)))) < 60]
#auc <- auc[apply(auc,1, function(x) length(which(is.na(x)))) < 20,]
auc[] <- auc %>% mutate_all(function(x) impute(x))                       # only for the survival prediction !


mut_original <- read.csv("dream_data_leaderboard/dnaseq.csv",sep=",") #CHECK
mut <- read.csv("features_validation/mut/dnaseq_full.csv")
mut$value <- 1
mut <- pivot_wider(data = mut[,c("lab_id","var_name","value")], names_from = var_name, values_from = value, names_repair = 'unique')
#mut <- spread(data = mut[,c("lab_id","var_name","value")], key = var_name, value = value)
mut <- as.data.frame(mut)
rownames <- mut$lab_id
mut[is.null(mut)] <- 0
mut[is.na(mut)] <- 0
row.names(mut) <- make.names(mut$lab_id)
mut <- mut[,-1]
mut <- mut[, apply(mut, 2, sum) >= 4 ]

samples <- row.names(rna) # bugfix
missing_mut <- samples[!samples %in% row.names(mut)]
missing_zeros <- matrix(0, nrow = length(missing_mut), ncol = ncol(mut))
row.names(missing_zeros) = missing_mut
mut <- rbind(mut, missing_zeros)

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
clin <- clin[,-ncol(clin)]
clin <- model.matrix(~., data=clin)

resp <- read.csv("features_validation/survival/response_full.csv",sep=",")
resp_original <- read.csv("dream_data_leaderboard/response.csv")
resp <- resp[,-1] # not in original
row.names(resp) = make.names(resp$lab_id)
resp <- resp[,-1]
resp$vitalStatus <- as.character(resp$vitalStatus)
resp$vitalStatus[resp$vitalStatus == "Alive"] <- 1
resp$vitalStatus[resp$vitalStatus == "Dead"] <- 0
resp$vitalStatus <- as.numeric(resp$vitalStatus)
resp <- resp[resp$overallSurvival != 0,]


#### Retrain models !
source("submission/sc1/input_data_functions.R")
auc <- import_aucs(path = "dream_data_leaderboard/aucs.csv")
rna <- import_rnaseq(path = "dream_data_leaderboard/rnaseq.csv")
mut <- import_dnaseq(path = "dream_data_leaderboard/dnaseq.csv")
clin <- import_clin(path_num = "dream_data_leaderboard/clinical_numerical.csv", path_cat = "dream_data_leaderboard/clinical_categorical.csv")



clin_feature_miss = setdiff(clin_feature, colnames(clin))
clin_feature_mat = matrix( data = 0, nrow = nrow(clin), ncol = length(clin_feature_miss),
                           dimnames = list(row.names(clin), clin_feature_miss))
clin = cbind(clin, clin_feature_mat)



directory <- "clin-surv"#"mut" #"rna"
descriptor <- directory
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
feature_path
response_path

save(clin, file = feature_path)
save(resp, file = response_path)





### DNN h2o offline prediction, should also work otherwise !

model <- h2o.loadModel(final_model_list$param[[1]][[1]][[1]])
modelfile <- h2o.download_mojo(model, path="metadata/h2odnn_pickle/", get_genmodel_jar=TRUE, genmodel_name = model@model_id)
object <- as.h2o(rna[,final_model_list$gene_names_filtered[[1]]%>%unlist])
h2o.saveModelDetails(object, path = "/storage/groups/cbm01/workspace/dream_aml/metadata/h2odnn_pickle/DeepLearning_model_R_1583746336351_1.json", force = F)


h2o.predict_json(model = "/storage/groups/cbm01/workspace/dream_aml/metadata/h2odnn_pickle/DeepLearning_model_R_1583746336351_1.zip",
                 genmodelpath = "/storage/groups/cbm01/workspace/dream_aml/metadata/h2odnn_pickle/DeepLearning_model_R_1583746336351_1",
                 )




