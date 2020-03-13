library(tidyr)
library(dplyr)
library(randomForestSRC)
library(randomForest)
library(survival)
library(Hmisc)

source("submission/sc2/input_data_functions.R")
source("submission/sc2/general.R")

mod_rna  <- loadRData("submission/sc2/models/rna-surv/rfsurv_default.RData")
mod_mut  <- loadRData("submission/sc2/models/mut-surv/rfsurv_default.RData")
mod_clin  <- loadRData("submission/sc2/models/clin-surv/rfsurv_default.RData")
mod_auc  <- loadRData("submission/sc2/models/auc-surv/rfsurv_default.RData")

rna <- import_rnaseq("dream_data/rnaseq.csv")
mut <- import_dnaseq("dream_data/dnaseq.csv")
clin <- import_clin(path_num = "dream_data/clinical_numerical.csv", 
                    path_cat = "dream_data/clinical_categorical.csv")
auc <- import_aucs("dream_data/aucs.csv")
clin_feature = mod_clin_glm$gene_names_filtered[[1]] %>% unlist()
clin_feature_miss = setdiff(clin_feature, colnames(clin))
clin_feature_mat = matrix( data = 0, nrow = nrow(clin), ncol = length(clin_feature_miss),
                           dimnames = list(row.names(clin), clin_feature_miss))
clin = cbind(clin, clin_feature_mat)

stacked_models <- loadRData("submission/sc2/models/stacked_models_sc2.RData")

id = intersect(rownames(auc),intersect(rownames(clin), intersect(rownames(mut), rownames(rna))))
clin = clin[id,]
mut = mut[id,]
rna = rna[id,]
auc = auc[id,]
  
predict_rna <- predict.rfsrc(mod_rna$param[[1]][[1]][[1]], newdata = as.data.frame(rna))$predicted; names(predict_rna)=row.names(rna)
predict_mut <- predict.rfsrc(mod_mut$param[[1]][[1]][[1]], newdata = as.data.frame(mut))$predicted; names(predict_mut)=row.names(mut)
predict_clin <- predict.rfsrc(mod_clin$param[[1]][[1]][[1]], newdata = as.data.frame(clin))$predicted; names(predict_clin)=row.names(clin)
predict_auc <- predict.rfsrc(mod_auc$param[[1]][[1]][[1]], newdata = as.data.frame(auc))$predicted; names(predict_auc)=row.names(auc)

pred_df = data.frame(mut.rf = predict_mut, rna.rf = predict_rna, clin.rf = predict_clin, auc.rf = predict_auc)
colnames(pred_df) = c("mut.rf","rna.rf","clin.rf","auc.rf")						
  
final_predict = predict.rfsrc(stacked_models$param[[1]][[1]][[1]], newdata = pred_df)$predicted; names(final_predict) = row.names(pred_df)


final_predict_df = data.frame(lab_id = names(final_predict), survival = final_predict)

final_predict_df$lab_id <- unlist(lapply(final_predict_df$lab_id, function(x) gsub("X","",x)))
final_predict_df$lab_id <- unlist(lapply(final_predict_df$lab_id, function(x) gsub("\\.","-",x)))
final_predict_df <- final_predict_df[,c("lab_id","survival")]
  

write.csv(final_predict_df, "predictions_sc2_test.csv", row.names = F)
