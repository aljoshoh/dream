#!/usr/bin/env Rscript
setwd("/storage/groups/cbm01/workspace/dream_aml/")
args <- as.numeric(commandArgs(trailingOnly = TRUE)); args <- args # get the submission script argument
print("Arguments:")
print(args)

## Alex' Models
## note that our docker image is using r3.6.1
source("R/learning.R")
source("R/run_pipeline.R")
source("R/feature_selection.R")
source("R/algorithms.R")
source("R/general.R")

#### SUBMISSION
# ./SUBMIT.sh alex_pipeline.R 8 5 <- old command
# ./SUBMIT.sh alex_pipeline.R 40 1
###############
#args <- (args[2]-1)*8+args[1] #8=number of jobs per array
args <- args[1] #<- new command 
print(paste0("Running with argument: ",as.character(args)))
numberofargs <- 8*5 # if sequential, set to 1
###############

directory <- "rna"#"mut" #"rna"
descriptor <- "rna"#"mut" #"rna"

feature_path = paste0("features_validation/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features_validation/",directory,"/",descriptor,"_response.RData") # path to response
rna <- get_features(feature_path)
auc <- get_features(response_path)
descriptor <- "glm"
#dump_features(rna, path = paste0("features/",directory,"/",descriptor,"_features.RData"))
#dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response.RData"))

if(numberofargs ==1){
  auc <- auc
} else {
  auc <- cut_df(auc, numberofargs,args)
}

dump_features(auc, path = paste0("features_validation/",directory,"/",descriptor,"_response_",as.character(args),".RData"))

models_list <- run_pipeline_benchmark(
  feature_path = paste0("features_validation/",directory,"/",descriptor,"_features.RData"), # path to features
  response_path = paste0("features_validation/",directory,"/",descriptor,"_response_",as.character(args),".RData"), # path to response
  submission = T,
  kfold = 10, 
  method = c(descriptor),
  hyperparam = c("alpha"=1.), # list(c(NULL),c(NULL)), #list(c(333),c(500)), # c("alpha"=1.),
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = 1,
  args = args,
  FUN = AnvSigNumFeature
  
)
# also possible to add FUN=AnvSigGen 
# @phong: the method "make_fit" does not yet return the results of the filtering

save(models_list, file = paste0("outputs/",directory,"/",descriptor,"_default._10fold_cvseed1_instance",as.character(args),".RData")) # "_test"





########### FURTHER STEPS FOR BUILDING THE FULL MODEL
# You can do these things also in a separate scripts, since these are not that computationally heavy
if(FALSE){

  cv_models <- import_cv_results(
    partial_path = paste0("dnn_","default._"),
    directory = paste0("outputs/",directory)
  )  
  
  
  lambda_min <- lambda_min(
    pipeline_object = models_list # the object created from the pipeline object, only works if returnFit=F
  )
  
  final_model_list <- run_pipeline_final(
    feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
    response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
    submission = T,
    FUN = function(x){return(x)},
    method = "dnn",
    hyperparam = NULL#c("alpha"=0.5)
  )
  
  save(final_models_list, file = paste0("outputs/",directory,"/","dnn_","default.RData"))
  
  
  # you can "predict(final_model_list[[1]], s = lambda_min[[1]], newx=blablala)" for choosing lambda with optimal cv-score
}
#####################################################






