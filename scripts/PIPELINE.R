#!/usr/bin/env Rscript
setwd("/storage/groups/cbm01/workspace/dream_aml/")
args <- as.numeric(commandArgs(trailingOnly = TRUE)); args <- args # get the submission script argument
print("Arguments:")
print(args)
#args <- 1 # TEMPORALY LOCALLY NEEDED FOR TESTING

## Alex' Models
## note that our docker image is using r3.6.1
source("R/learning.R")
source("R/run_pipeline.R")
source("R/feature_selection.R")
source("R/algorithms.R")
source("R/general.R")

#### how to run on cluster
# ./SUBMIT.sh PIPELINE.R numberofargs 1
args <- args[1] #<- new command 
print(paste0("Running with argument: ",as.character(args)))
##########################

### SCRIPT PARAMETER
directory <- "clin-auc"#"mut" #"rna"
descriptor <- "dnn" # the descriptor means the method in this script, not the same as in PREPROCESS.R
param <- list(c(NULL),c(NULL)) #list(c(333),c(500)) # c("alpha"=1.)
####################
if(descriptor=="dnn"){h2o.init(port=8506)}
models_list <- run_pipeline_benchmark(
  feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
  response_path = paste0("features/",directory,"/",descriptor,"_response_",as.character(args),".RData"), # path to response
  submission = T,
  kfold = 10, 
  method = c(descriptor),
  hyperparam = param,
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = 1,
  FUN = AnvSigNumFeature,#AnvSigSurvFeature,
  args = args
)
save(models_list, file = paste0("outputs/",directory,"/",descriptor,"_default._10fold_cvseed1_instance",as.character(args),".RData"))


h2o.shutdown()


########### FURTHER STEPS FOR BUILDING THE FULL MODEL
# You can do these things also in a separate scripts, since these are not that computationally heavy
if(FALSE){

  cv_models <- import_cv_results(
    partial_path = paste0(descriptor,"_default._"),
    directory = paste0("outputs/",directory)
  )  
  
  save(cv_models, file = paste0("outputs/",directory,"/",descriptor,"_default_cv.RData"))
  
  if(descriptor=="glm"){
    lambda.min <- lambda_min(
      pipeline_object = cv_models # the object created from the pipeline object, only works if returnFit=F
    )
    save(lambda.min, file =  paste0("outputs/",directory,"/",descriptor,"_default_lambda.min.RData"))
  }
  
  if(descriptor=="dnn"){h2o.init(port=8510)}
  final_model_list <- run_pipeline_final(
    feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
    response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
    submission = T,
    method = descriptor,
    hyperparam = param,
    FUN = AnvSigNumFeature #AnvSigSurvFeature
  )
  
  save(final_model_list, file = paste0("outputs/",directory,"/",descriptor,"_default.RData"))
  
  
  # you can "predict(final_model_list[[1]], s = lambda_min[[1]], newx=blablala)" for choosing lambda with optimal cv-score
}
#####################################################






