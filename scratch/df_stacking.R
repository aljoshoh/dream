## Random forest model stacking


## 
library(randomForest)
library(caret)
source("R/learning.R")
source("R/run_pipeline.R")
source("R/select_gene_sc1.R")
source("R/algorithms.R")


feature_path = "features/alex_features.RData" # path to features
response_path = "features/alex_phenotypes.RData" # path to response
rna <- get_features(feature_path)
auc <- get_features(response_path)

rna <- rna[,1:1000]
auc <- auc[,1:2]

dump_features(auc, path = "features/alex_phenotypes_red.RData")
dump_features(rna, path = "features/alex_features_red.RData")
dump_features(models_list, path = "features/alex_rf_models.RData")


modelsa <- get_features("features/alex_glm_models.RData", matrixfy = F)
modelsb <- get_features("features/alex_rf_models.RData", matrixfy = F)


### $cv of the two benchmarking objects must match

preds_stack <- lapply(1:ncol(auc), function(y){
model1 <- unlist(lapply(1:length(modelsb$cv$test_sets), function(x) predict(modelsb$param[[x,y]][[1]], newdata = rna[modelsb$cv$test_sets[[x]],])));
model2 <- unlist(lapply(1:length(modelsa$cv$test_sets), function(x) predict(modelsa$param[[x,y]][[1]], newx = rna[modelsa$cv$test_sets[[x]],], s = 'lambda.min')[,1]));
temp <- cbind(model1, model2);
return(temp)
})

stack_features <- preds_stack
dump_features(stack_features, path = "features/alex_stacked_models_test.RData")

models_list_stacked <- run_pipeline_benchmark(
  feature_path = "features/alex_stacked_models_test.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
  response_path = "features/alex_phenotypes_red.RData", # path to response
  submission = F,
  kfold = NULL, 
  method = "glm",
  hyperparam = c("alpha"=0.5), #list(c(333),c(500)), # c("alpha"=0.5),
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = NULL, # supply the parallel processing counter
  CVBuilt = modelsa$cv,
  stack = T
)
models_list_stacked$score












