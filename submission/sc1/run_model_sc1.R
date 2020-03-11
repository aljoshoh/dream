library(tidyr)
library(dplyr)
library(randomForest)
library(survival)
library(glmnet)
library(h2o)
library(readr)
library(Hmisc)

source("/home/input_data_functions.R") ### alex bug fix, was no absolute path
source("/home/general.R")

rna <- import_rnaseq("/input/rnaseq.csv")
mut <- import_dnaseq("/input/dnaseq.csv")
clin <- import_clin(path_num = "/input/clinical_numerical.csv", path_cat = "/input/clinical_categorical.csv") 
 
mod_rna_glm  <- loadRData("/home/models/rna-auc/glm_default.RData")
mod_rna_dnn  <- loadRData("/home/models/rna-auc/dnn_default.RData")
mod_rna_rf  <- loadRData("/home/models/rna-auc/rf_default.RData")
mod_mut_glm  <- loadRData("/home/models/mut-auc/glm_default.RData")
mod_mut_dnn  <- loadRData("/home/models/mut-auc/dnn_default.RData")
mod_mut_rf  <- loadRData("/home/models/mut-auc/rf_default.RData")
mod_clin_glm  <- loadRData("/home/models/clin-auc/glm_default.RData")
mod_clin_dnn  <- loadRData("/home/models/clin-auc/dnn_default.RData")
mod_clin_rf  <- loadRData("/home/models/clin-auc/rf_default.RData")

inhibitor_name <- loadRData("/home/drug_names.RData")

lambda.min_rna <- loadRData("/home/models/rna-auc/glm_default_lambda.min.RData")
lambda.min_mut <- loadRData("/home/models/mut-auc/glm_default_lambda.min.RData")
lambda.min_clin <- loadRData("/home/models/clin-auc/glm_default_lambda.min.RData")

h2o.init()
PATH = "/home/models/h2o_models/"
id = intersect(rownames(clin), intersect(rownames(mut), rownames(rna)))
clin = clin[id,]
mut = mut[id,]
rna = rna[id,]

final_predict = list()

for (i in 1:length(inhibitor_name)){

  lambda_rna = lambda.min_rna[[i]]
  lambda_mut = lambda.min_mut[[i]]
  lambda_clin = lambda.min_clin[[i]]

  mod_mut_dnn$param[[i]][[1]][[1]] = convert_path(path = PATH, path.in.object = mod_mut_dnn$param[[i]][[1]][[1]])
  mod_rna_dnn$param[[i]][[1]][[1]] = convert_path(path = PATH, path.in.object = mod_rna_dnn$param[[i]][[1]][[1]])
  mod_clin_dnn$param[[i]][[1]][[1]] = convert_path(path = PATH, path.in.object = mod_clin_dnn$param[[i]][[1]][[1]])

  predict_rna_glm = predict(mod_rna_glm$param[[i]][[1]][[1]], rna[,mod_rna_glm$gene_names_filtered[[i]] %>% unlist], s = lambda_rna)
  predict_rna_dnn = predict(h2o.loadModel(mod_rna_dnn$param[[i]][[1]][[1]]), newdata = as.h2o(rna[,mod_rna_dnn$gene_names_filtered[[i]] %>% unlist])) %>% as.data.frame%>%unlist
  names(predict_rna_dnn) = rownames(rna)
  predict_rna_rf = predict(mod_rna_rf$param[[i]][[1]][[1]], rna[,mod_rna_rf$gene_names_filtered[[i]] %>% unlist])
  predict_mut_glm = predict(mod_mut_glm$param[[i]][[1]][[1]], mut[,mod_mut_glm$gene_names_filtered[[i]] %>% unlist], s = lambda_mut)
  predict_mut_dnn = predict(h2o.loadModel(mod_mut_dnn$param[[i]][[1]][[1]]), newdata = as.h2o(mut[,mod_mut_dnn$gene_names_filtered[[i]] %>% unlist])) %>% as.data.frame%>%unlist
  names(predict_mut_dnn) = rownames(mut)
  predict_mut_rf = predict(mod_mut_rf$param[[i]][[1]][[1]], mut[,mod_mut_rf$gene_names_filtered[[i]] %>% unlist])
  predict_clin_glm = predict(mod_clin_glm$param[[i]][[1]][[1]], clin[,mod_clin_glm$gene_names_filtered[[i]] %>% unlist], s = lambda_clin)
  predict_clin_dnn = predict(h2o.loadModel(mod_clin_dnn$param[[i]][[1]][[1]]), newdata = as.h2o(clin[,mod_clin_glm$gene_names_filtered[[i]] %>% unlist])) %>% as.data.frame%>%unlist
  names(predict_clin_dnn) = rownames(clin)
  predict_clin_rf = predict(mod_clin_rf$param[[i]][[1]][[1]], clin[,mod_clin_rf$gene_names_filtered[[i]] %>% unlist])

  pred_df = data.frame(mut.glm = predict_mut_glm, mut.rf = predict_mut_rf, mut.dnn = predict_mut_dnn,
						rna.glm = predict_rna_glm, rna.rf = predict_rna_rf, rna.dnn = predict_rna_dnn,
						clin.glm = predict_clin_glm, clin.rf = predict_clin_rf, clin.dnn = predict_clin_dnn, row.names = names(predict_rna_rf))
  colnames(pred_df) = c("mut.glm","mut.rf","mut.dnn","rna.glm","rna.rf","rna.dnn","clin.glm","clin.rf","clin.dnn")						

  final_predict[[i]] = predict(stacked_models$param[[i]][[1]][[1]], newdata = pred_df)

}

names(final_predict) = inhibitor_name
final_predict_vec = unlist(final_predict)

final_predict_df = data.frame(identifier = names(final_predict_vec), auc = final_predict_vec,
                              row.names = NULL) %>%
  separate(col = identifier, into = c("inhibitor", "lab_id1", "labid2"), sep = "\\.") %>%
  unite(col = "labid", 2:3, sep = ".")
  

write.csv(final_predict_df, "/output/predictions.csv", row.names = F)

