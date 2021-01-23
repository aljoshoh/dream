models1 <- get_features("submission/sc1/models/mut-auc/dnn_default.RData", matrixfy = F)
models2 <- get_features("submission/sc1/models/rna-auc/dnn_default.RData", matrixfy = F)
models3 <- get_features("submission/sc1/models/clin-auc/dnn_default.RData", matrixfy = F)

dnn_models_sc1 <- c(models1$param %>% unlist,
  models2$param %>% unlist,
  models3$param %>% unlist)

for(model in dnn_models_sc1){
  print(model)
  system(paste0("cp ",model," /storage/groups/cbm01/workspace/dream_aml/submission/sc1/models/h2o_models"))
}
