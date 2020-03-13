#Thanks Phong for crafting this !:)
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(randomForest)
  library(survival)
  library(glmnet)
  library(h2o)
  library(readr)
  library(Hmisc)
  library(caret)
  library(glmnet)
})

source("/usr/local/bin/input_data_functions.R")
source("/usr/local/bin/general.R")

rna <- import_rnaseq("/input/rnaseq.csv")
mut <- import_dnaseq("/input/dnaseq.csv")
clin <- import_clin(path_num = "/input/clinical_numerical.csv", 
                    path_cat = "/input/clinical_categorical.csv")
auc <- import_aucs("/input/aucs.csv")

mod_rna  <- loadRData("/usr/local/bin/models/rna-surv/rfsurv_default.RData")
mod_mut  <- loadRData("/usr/local/bin/models/mut-surv/rfsurv_default.RData")   # alex edits: changing /home to /usr/local/bin
mod_clin  <- loadRData("/usr/local/bin/models/clin-surv/rfsurv_default.RData")
mod_auc  <- loadRData("/usr/local/bin/models/auc-surv/rfsurv_default.RData")

stacked_models <- loadRData("/usr/local/bin/models/stacked_models_sc2.RData")

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
  

write.csv(final_predict_df, "/output/predictions.csv", row.names = F)

