## Random forest model stacking


## 
library(randomForest)
library(caret)
source("R/learning.R")
source("R/run_pipeline.R")
source("R/select_gene_sc1.R")
source("R/algorithms.R")
source("R/general.R")

#feature_path = "features/alex_features.RData" # path to features
#response_path = "features/alex_phenotypes.RData" # path to response
#rna <- get_features(feature_path)
#auc <- get_features(response_path)

#rna <- rna[,1:1000]
#auc <- auc[,1:2]

#dump_features(auc, path = "features/alex_phenotypes_red.RData")
#dump_features(rna, path = "features/alex_features_red.RData")
#dump_features(models_list, path = "features/alex_rf_models.RData")


modelsa <- get_features("outputs/mut-auc/glm_default_cv.RData", matrixfy = F)
modelsb <- get_features("outputs/mut-auc/rf_default_cv.RData", matrixfy = F)
modelsc <- get_features("outputs/mut-auc/dnn_default_cv.RData", matrixfy = F)
modelsd <- get_features("outputs/rna-auc/glm_default_cv.RData", matrixfy = F)
#modelse <- get_features("outputs/rna-auc/rf_default_cv.RData", matrixfy = F)
#modelsf <- get_features("outputs/rna-auc/dnn_default_cv.RData", matrixfy = F)
modelsg <- get_features("outputs/clin-auc/glm_default_cv.RData", matrixfy = F)
modelsh <- get_features("outputs/clin-auc/rf_default_cv.RData", matrixfy = F)
modelsi <- get_features("outputs/clin-auc/dnn_default_cv.RData", matrixfy = F)

h2o.init()
modelsc$param_paths <- modelsc$param
#modelsf$param_paths <- modelsc$param
modelsi$param_paths <- modelsc$param

### SCRIPT PARAMETER
directory <- "rna-auc"#"mut" #"rna"
descriptor <- "glm" # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
rna <- get_features(feature_path)
auc_rna <- get_features(response_path)

### SCRIPT PARAMETER
directory <- "mut-auc"#"mut" #"rna"
descriptor <- "glm" # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
mut <- get_features(feature_path)
auc_mut <- get_features(response_path)

### SCRIPT PARAMETER
directory <- "clin-auc"#"mut" #"rna"
descriptor <- "glm" # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
clin <- get_features(feature_path)
auc_clin <- get_features(response_path)

#DEPRECATED
#preds_stack <- lapply(1:ncol(auc), function(y){
#model1 <- unlist(lapply(1:length(modelsa$cv$test_sets), function(x) predict(modelsa$param[[x,y]][[1]], newx = rna[modelsa$cv$test_sets[[x]],], s = 'lambda.min')[,1]));
#model2 <- unlist(lapply(1:length(modelsb$cv$test_sets), function(x) predict(modelsb$param[[x,y]][[1]], newdata = rna[modelsb$cv$test_sets[[x]],])));
#model3 <- unlist(lapply(1:length(modelsc$cv$test_sets), function(x) predict(h2o.loadModel(modelsc$param[[x,y]][[1]]), newdata = as.h2o(rna[modelsc$cv$test_sets[[x]],]))%>%as.data.frame%>%unlist));
#temp <- cbind(model1, model2, model3);
#return(temp)
#})


preds_stack <- list()
for(y in 1:ncol(auc)){  #1:ncol(auc)
  model1 <- unlist(lapply(1:length(modelsa$cv$test_sets), function(x) predict(modelsa$param[[x,y]][[1]], newx = mut[modelsa$cv$test_sets[[x]],], s = 'lambda.min')[,1]))
  model2 <- unlist(lapply(1:length(modelsb$cv$test_sets), function(x) predict(modelsb$param[[x,y]][[1]], newdata = mut[modelsb$cv$test_sets[[x]],])))
  model3 <- unlist(lapply(1:length(modelsc$cv$test_sets), function(x) predict(h2o.loadModel(modelsc$param[[x,y]][[1]]), newdata = as.h2o(mut[modelsc$cv$test_sets[[x]],]))%>%as.data.frame%>%unlist))
  model4 <- unlist(lapply(1:length(modelsd$cv$test_sets), function(x) predict(modelsd$param[[x,y]][[1]], newx = rna[modelsd$cv$test_sets[[x]],modelsd$gene_names_filtered[[x,y]][[1]]], s = 'lambda.min')[,1]))
  #model5 <- unlist(lapply(1:length(modelse$cv$test_sets), function(x) predict(modelse$param[[x,y]][[1]], newdata = rna[modelse$cv$test_sets[[x]],modele$gene_names_filtered[[x,y]][[1]]])))
  #model6 <- unlist(lapply(1:length(modelsf$cv$test_sets), function(x) predict(h2o.loadModel(modelsf$param[[x,y]][[1]]), newdata = as.h2o(rna[modelsf$cv$test_sets[[x]],modelsd$gene_names_filtered[[x,y]][[1]]]))%>%as.data.frame%>%unlist))
  model7 <- unlist(lapply(1:length(modelsg$cv$test_sets), function(x) predict(modelsg$param[[x,y]][[1]], newx = clin[modelsg$cv$test_sets[[x]],], s = 'lambda.min')[,1]))
  model8 <- unlist(lapply(1:length(modelsh$cv$test_sets), function(x) predict(modelsh$param[[x,y]][[1]], newdata = clin[modelsh$cv$test_sets[[x]],])))
  model9 <- unlist(lapply(1:length(modelsi$cv$test_sets), function(x) predict(h2o.loadModel(modelsi$param[[x,y]][[1]]), newdata = as.h2o(clin[modelsi$cv$test_sets[[x]],]))%>%as.data.frame%>%unlist))
  intsec <- intersect(intersect(names(model1),names(model4)), names(model7))
  names(model3) = names(model1)
  #names(model6) = names(model4)
  names(model9) = names(model7)
  temp <- cbind(model1[intsec], 
                model2[intsec], 
                model3[intsec],
                model4[intsec],
                #model5[intsec],
                #model6[intsec],
                model7[intsec],
                model8[intsec],
                model9[intsec]
                );
  preds_stack[[y]] <- temp
  print(paste0(as.character(y/ncol(auc)*100),"%"))
}

stack_features <- preds_stack
dump_features(stack_features, path = "features/stacked_models_sc1.RData")
models_list_stacked <- run_pipeline_benchmark(
  feature_path = "features/stacked_models_sc1.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
  response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
  submission = F,
  kfold = NULL, 
  method = "rf",
  hyperparam = list(c(NULL),c(NULL)), #c("alpha"=0.5), #list(c(333),c(500)), # c("alpha"=0.5),
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = NULL, # supply the parallel processing counter
  CVBuilt = modelsa$cv,
  stack = T
)

models_list_stacked <- run_pipeline_final(
  feature_path = "features/stacked_models_sc1.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
  response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
  submission = T,
  method = "rf",
  hyperparam = list(c(NULL),c(NULL)), #c("alpha"=0.5), #list(c(333),c(500)), # c("alpha"=0.5),
  stack = T
)

dump_features(models_list_stacked, path = "outputs/stacked_models_sc1.RData")




### FUNCTIONS
function(model, x, y){
  model <- h2o.loadModel(model$param[[x,y]][[1]])
  return(out)
}




