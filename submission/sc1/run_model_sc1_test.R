lib = "/storage/groups/cbm01/workspace/phong.nguyen_new/Rlib/"

library(tidyr)
library(dplyr)
library(randomForest, lib.loc = lib)
library(survival, lib.loc = lib)
library(glmnet, lib.loc = lib)
library(h2o, lib.loc = lib)
library(Hmisc, lib.loc = lib)

source("submission/sc1/input_data_functions.R")
source("R/general.R")


mod_rna_glm  <- loadRData("submission/sc1/models/rna-auc/glm_default.RData")
mod_rna_dnn  <- loadRData("submission/sc1/models/rna-auc/dnn_default.RData")
mod_rna_rf  <- loadRData("submission/sc1/models/rna-auc/rf_default.RData")
mod_mut_glm  <- loadRData("submission/sc1/models/mut-auc/glm_default.RData")
mod_mut_dnn  <- loadRData("submission/sc1/models/mut-auc/dnn_default.RData")
mod_mut_rf  <- loadRData("submission/sc1/models/mut-auc/rf_default.RData")
mod_clin_glm  <- loadRData("submission/sc1/models/clin-auc/glm_default.RData")
mod_clin_dnn  <- loadRData("submission/sc1/models/clin-auc/dnn_default.RData")
mod_clin_rf  <- loadRData("submission/sc1/models/clin-auc/rf_default.RData")


rna <- import_rnaseq("dream_data/rnaseq.csv")
mut <- import_dnaseq("dream_data_leaderboard/dnaseq.csv")
clin <- import_clin(path_num = "dream_data_leaderboard/clinical_numerical.csv", 
                    path_cat = "dream_data_leaderboard/clinical_categorical.csv")
clin_feature = mod_clin_glm$gene_names_filtered[[1]] %>% unlist()
clin_feature_miss = setdiff(clin_feature, colnames(clin))
clin_feature_mat = matrix( data = 0, nrow = nrow(clin), ncol = length(clin_feature_miss),
                           dimnames = list(row.names(clin), clin_feature_miss))
clin = cbind(clin, clin_feature_mat)


stacked_models <- loadRData("submission/sc1/models/stacked_models_sc1.RData")

inhibitor_name <- loadRData("submission/sc1/drug_names.RData")

lambda.min_rna <- loadRData("submission/sc1/models/rna-auc/glm_default_lambda.min.RData")
lambda.min_mut <- loadRData("submission/sc1/models/mut-auc/glm_default_lambda.min.RData")
lambda.min_clin <- loadRData("submission/sc1/models/clin-auc/glm_default_lambda.min.RData")

final_predict = list()

h2o.init()
PATH = "/storage/groups/cbm01/workspace/dream_aml/submission/sc1/models/h2o_models/"
id = intersect(rownames(clin), intersect(rownames(mut), rownames(rna)))
clin = clin[id,]
mut = mut[id,]
rna = rna[id,]

for (i in 1:length(inhibitor_name)) {
  
  lambda_rna = lambda.min_rna[[i]]
  lambda_mut = lambda.min_mut[[i]]
  lambda_clin = lambda.min_clin[[i]]
  
  mod_mut_dnn$param[[i]][[1]][[1]] = convert_path(path = PATH, path.in.object = mod_mut_dnn$param[[i]][[1]][[1]])
  mod_rna_dnn$param[[i]][[1]][[1]] = convert_path(path = PATH, path.in.object = mod_rna_dnn$param[[i]][[1]][[1]])
  mod_clin_dnn$param[[i]][[1]][[1]] = convert_path(path = PATH, path.in.object = mod_clin_dnn$param[[i]][[1]][[1]])
  
  predict_rna_glm = predict(mod_rna_glm$param[[i]][[1]][[1]], rna[,mod_rna_glm$gene_names_filtered[[i]] %>% unlist], s = lambda_rna)
  predict_rna_dnn = predict(h2o.loadModel(mod_rna_dnn$param[[i]][[1]][[1]]), 
                            newdata = as.h2o(rna[,mod_rna_dnn$gene_names_filtered[[i]] %>% unlist])) %>% as.data.frame%>%unlist
  names(predict_rna_dnn) = rownames(rna)
  predict_rna_rf = predict(mod_rna_rf$param[[i]][[1]][[1]], rna[,mod_rna_rf$gene_names_filtered[[i]] %>% unlist])
  predict_mut_glm = predict(mod_mut_glm$param[[i]][[1]][[1]], mut[,mod_mut_glm$gene_names_filtered[[i]] %>% unlist], s = lambda_mut)
  predict_mut_dnn = predict(h2o.loadModel(mod_mut_dnn$param[[i]][[1]][[1]]), 
                            newdata = as.h2o(mut[,mod_mut_dnn$gene_names_filtered[[i]] %>% unlist])) %>% as.data.frame%>%unlist
  names(predict_mut_dnn) = rownames(mut)
  predict_mut_rf = predict(mod_mut_rf$param[[i]][[1]][[1]], mut[,mod_mut_rf$gene_names_filtered[[i]] %>% unlist])
  predict_clin_glm = predict(mod_clin_glm$param[[i]][[1]][[1]], clin[,mod_clin_glm$gene_names_filtered[[i]] %>% unlist], s = lambda_clin)
  predict_clin_dnn = predict(h2o.loadModel(mod_clin_dnn$param[[i]][[1]][[1]]), 
                             newdata = as.h2o(clin[,mod_clin_glm$gene_names_filtered[[i]] %>% unlist])) %>% as.data.frame%>%unlist
  names(predict_clin_dnn) = rownames(clin)
  predict_clin_rf = predict(mod_clin_rf$param[[i]][[1]][[1]], clin[,mod_clin_rf$gene_names_filtered[[i]] %>% unlist])
  
  pred_df = data.frame(mut.glm = predict_mut_glm, mut.rf = predict_mut_rf, mut.dnn = predict_mut_dnn,
                       rna.glm = predict_rna_glm, rna.rf = predict_rna_rf, rna.dnn = predict_rna_dnn,
                       clin.glm = predict_clin_glm, clin.rf = predict_clin_rf, clin.dnn = predict_clin_dnn)
  colnames(pred_df) = c("mut.glm","mut.rf","mut.dnn","rna.glm","rna.rf","rna.dnn",
                        "clin.glm","clin.rf","clin.dnn")						
  
  final_predict[[i]] = predict(stacked_models$param[[i]][[1]][[1]], newdata = pred_df)
}

names(final_predict) = inhibitor_name
final_predict_vec = unlist(final_predict)

final_predict_df = data.frame(identifier = names(final_predict_vec), auc = final_predict_vec,
                              row.names = NULL) %>%
  separate(col = identifier, into = c("inhibitor", "lab_id1", "labid2"), sep = "\\.") %>%
  unite(col = "labid", 2:3, sep = ".")
  

write.csv(final_predict_df, "predictions_test.csv", row.names = F)
