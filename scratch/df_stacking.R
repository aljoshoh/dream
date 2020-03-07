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
#modelsd <- get_features("outputs/rna-auc/glm_default_cv.RData", matrixfy = F)
#modelse <- get_features("outputs/rna-auc/rf_default_cv.RData", matrixfy = F)
#modelsf <- get_features("outputs/rna-auc/dnn_default_cv.RData", matrixfy = F)

h2o.init()
modelsc$param_paths <- modelsc$param

### SCRIPT PARAMETER
directory <- "mut-auc"#"mut" #"rna"
descriptor <- "dnn" # the descriptor means the method in this script, not the same as in PREPROCESS.R
####################

### $cv of the two benchmarking objects must match
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
rna <- get_features(feature_path)
auc <- get_features(response_path)

#DEPRECATED
#preds_stack <- lapply(1:ncol(auc), function(y){
#model1 <- unlist(lapply(1:length(modelsa$cv$test_sets), function(x) predict(modelsa$param[[x,y]][[1]], newx = rna[modelsa$cv$test_sets[[x]],], s = 'lambda.min')[,1]));
#model2 <- unlist(lapply(1:length(modelsb$cv$test_sets), function(x) predict(modelsb$param[[x,y]][[1]], newdata = rna[modelsb$cv$test_sets[[x]],])));
#model3 <- unlist(lapply(1:length(modelsc$cv$test_sets), function(x) predict(h2o.loadModel(modelsc$param[[x,y]][[1]]), newdata = as.h2o(rna[modelsc$cv$test_sets[[x]],]))%>%as.data.frame%>%unlist));
#temp <- cbind(model1, model2, model3);
#return(temp)
#})

preds_stack <- list()
for(y in 1:ncol(auc)){
  model1 <- unlist(lapply(1:length(modelsa$cv$test_sets), function(x) predict(modelsa$param[[x,y]][[1]], newx = rna[modelsa$cv$test_sets[[x]],], s = 'lambda.min')[,1]))
  model2 <- unlist(lapply(1:length(modelsb$cv$test_sets), function(x) predict(modelsb$param[[x,y]][[1]], newdata = rna[modelsb$cv$test_sets[[x]],])))
  model3 <- unlist(lapply(1:length(modelsc$cv$test_sets), function(x) predict(h2o.loadModel(modelsc$param[[x,y]][[1]]), newdata = as.h2o(rna[modelsc$cv$test_sets[[x]],]))%>%as.data.frame%>%unlist))
  temp <- cbind(model1, model2, model3);
  preds_stack[[y]] <- temp
  print(paste0(as.character(y/ncol(auc)*100),"%"))
}

stack_features <- preds_stack
dump_features(stack_features, path = "features/alex_stacked_models_muttest.RData")
models_list_stacked <- run_pipeline_benchmark(
  feature_path = "features/alex_stacked_models_muttest.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
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
  feature_path = "features/alex_stacked_models_muttest.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
  response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
  submission = T,
  method = "rf",
  hyperparam = list(c(NULL),c(NULL)), #c("alpha"=0.5), #list(c(333),c(500)), # c("alpha"=0.5),
  stack = T
)

dump_features(models_list_stacked, path = "outputss/alex_stacked_models_muttest.RData")




### FUNCTIONS
function(model, x, y){
  model <- h2o.loadModel(model$param[[x,y]][[1]])
  return(out)
}




