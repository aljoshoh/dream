### SC1
library(glmnet)
library(h2o)
library(dplyr)
### SCRIPT PARAMETER
directory <- ""#"mut-surv" #"rna-auc"
descriptor <- "stacked_models_sc1"#"cox"#"stacked_models_sc1" # the descriptor does not care here #"rf"
model_object  <- loadRData(paste0("outputs/",directory,"/",descriptor,"_default_cv_classes.RData"))
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
feat <- get_features(feature_path)
resp <- get_features(response_path)

if(descriptor == "dnn"){h2o.init(ip = "localhost")}
preds_stack <- list()
for(y in 1:ncol(resp)){
  if(descriptor == "glm"){model <- unlist(lapply(1:length(model_object$cv$test_sets), function(x) predict(model_object$param[[x,y]][[1]], newx = feat[model_object$cv$test_sets[[x]],model_object$gene_names_filtered[[x,y]][[1]]], s = 'lambda.min')[,1]))}
  if(descriptor == "rf"){model <- unlist(lapply(1:length(model_object$cv$test_sets), function(x) predict(model_object$param[[x,y]][[1]], newdata = feat[model_object$cv$test_sets[[x]],model_object$gene_names_filtered[[x,y]][[1]]])))}
  if(descriptor == "dnn"){model <- unlist(lapply(1:length(model_object$cv$test_sets), function(x) predict(h2o.loadModel(model_object$param[[x,y]][[1]]), newdata = as.h2o(feat[model_object$cv$test_sets[[x]],model_object$gene_names_filtered[[x,y]][[1]]]))%>%as.data.frame%>%unlist))}
  if(descriptor == "stacked_models_sc1"){model <- unlist(lapply(1:length(model_object$cv$test_sets), function(x) predict(model_object$param[[x,y]][[1]], newdata = feat[[y]][model_object$cv$test_sets[[x]],model_object$gene_names_filtered[[x,y]][[1]]])))}
  if(descriptor == "cox"){model <- unlist(lapply(1:length(model_object$cv$test_sets), function(x) predict(model_object$param[[x,y]][[1]], newx = feat[model_object$cv$test_sets[[x]],model_object$gene_names_filtered[[x,y]][[1]]], s = 'lambda.min')[,1]))}
  preds_stack[[y]] <- model
  print(paste0(as.character(y/ncol(resp)*100),"%"))
}


performance <- unlist(lapply(1:ncol(resp), function(y){
  tmp <- mean(unlist(lapply(1:10, function(x){names(preds_stack[[y]])=unlist(model_object$cv$test_sets);cor(preds_stack[[y]][model_object$cv$test_sets[[x]]], resp[model_object$cv$test_sets[[x]],y], use = "complete.obs", method = "p")})));
  if(is.na(tmp)){tmp <- 0}else{tmp}
}))
mean(performance, na.rm = T)
save(performance, file = paste0("outputs/performances/",directory,"_",descriptor,"_performance.RData"))


if(FALSE){
  list <- list.files("/storage/groups/cbm01/workspace/dream_aml/outputs/performances/", pattern=".RData", all.files=FALSE, #meth, gex
             full.names=TRUE)
  collection <- c("stacked","clin.dnn","clin.glm","clin.rf","mut.dnn","mut.glm","mut.rf","rna.dnn","rna.glm","rna.rf")
  array <- matrix(NA, nrow = 122, ncol = 10)
  for(i in 1:length(list)){
    array[,i] <- loadRData(list[i])
  }
  lapply(1:122, function(x) max(na.omit(array[x,]))) %>% unlist %>% mean
  which.model <- lapply(1:122, function(x) collection[which.max(na.omit(array[x,]))]) %>% unlist
}
save(which.model, file = paste0("outputs/performances/",directory,"_",descriptor,"_which.model.RData"))

performance <- unlist(lapply(1:ncol(resp), function(y){
  tmp <- mean(unlist(lapply(1:10, function(x){model_object$score[[x,y]]})));
  if(is.na(tmp)){tmp <- 0}else{tmp}
}))
mean(performance, na.rm = T)
