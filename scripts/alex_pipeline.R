#!/usr/bin/env Rscript
setwd("/storage/groups/cbm01/workspace/dream_aml/")
args <- as.numeric(commandArgs(trailingOnly = TRUE)); args <- args # get the submission script argument
print("Arguments:")
print(args)

## Alex' Models
## note that our docker image is using r3.6.1
source("R/learning.R")
source("R/run_pipeline.R")
source("R/select_gene_sc1.R")
source("R/algorithms.R")
source("R/general.R")

#### SUBMISSION
# ./SUBMIT.sh alex_pipeline.R 8 5
###############
args <- (args[2]-1)*8+args[1] #8=number of jobs per array
print(paste0("Running with argument: ",as.character(args)))
numberofargs <- 8*5 # if sequential, set to 1
###############

feature_path = "features/mut/mut_features.RData" # path to features
response_path = "features/mut/mut_response.RData" # path to response
rna <- get_features(feature_path)
auc <- get_features(response_path)
dump_features(rna, path = "features/mut/mut_features.RData")
dump_features(auc, path = "features/mut/mut_response.RData")


auc <- cut_df(auc, numberofargs,args)	
dump_features(auc, path = paste0("features/mut/mut_response_",as.character(args),".RData"))

models_list <- run_pipeline_benchmark(
  feature_path = "features/mut/mut_features.RData", # path to features
  response_path = paste0("features/mut/mut_response_",as.character(args),".RData"), # path to response
  submission = T,
  kfold = 10, 
  method = c("dnn"),
  hyperparam = NULL, #list(c(NULL),c(NULL)), #list(c(333),c(500)), # c("alpha"=0.5),
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = 1 #args # supply the parallel processing counter
)
# also possible to add FUN=AnvSigGen 
# @phong: the method "make_fit" does not yet return the results of the filtering

save(models_list, file = paste0("outputs/mut/","dnn_","default._10fold_cvseed1_instance",as.character(args),".RData"))





########### FURTHER STEPS FOR BUILDING THE FULL MODEL
# You can do these things also in a separate scripts, since these are not that computationally heavy
if(FALSE){
  
  lambda_min <- lambda_min(
    pipeline_object = models_list # the object created from the pipeline object, only works if returnFit=F
  )
  
  final_model_list <- run_pipeline_final(
    feature_path = "features/alex_features.RData", # path to features
    response_path = "features/alex_phenotypes.RData", # path to response
    submission = T,
    FUN = function(x){return(x)},
    method = "glm",
    hyperparam = c("alpha"=0.5)
  )
  # you can "predict(final_model_list[[1]], s = lambda_min[[1]], newx=blablala)" for choosing lambda with optimal cv-score
}
#####################################################
