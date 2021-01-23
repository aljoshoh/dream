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

directory <- "mut"#"mut" #"rna"
descriptor <- "mut"#"mut" #"rna"

feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
rna <- get_features(feature_path)
auc <- get_features(response_path)
dump_features(rna, path = paste0("features/",directory,"/",descriptor,"_features.RData"))
dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response.RData"))


auc <- cut_df(auc, numberofargs,args)	
dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response_",as.character(args),".RData"))

models_list <- run_pipeline_benchmark(
  feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
  response_path = paste0("features/",directory,"/",descriptor,"_response_",as.character(args),".RData"), # path to response
  submission = T,
  kfold = 10, 
  method = c("glm"),
  hyperparam = c("alpha"=0.5), #list(c(NULL),c(NULL)), #list(c(333),c(500)), # c("alpha"=0.5),
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = 1,
  args = args
  #FUN = AnvSigNumFeature
  
)
# also possible to add FUN=AnvSigGen 
# @phong: the method "make_fit" does not yet return the results of the filtering

save(models_list, file = paste0("outputs/",directory,"/","dnn_","default._10fold_cvseed1_instance",as.character(args),"_test.RData"))





########### FURTHER STEPS FOR BUILDING THE FULL MODEL
# You can do these things also in a separate scripts, since these are not that computationally heavy
if(FALSE){

  cv_models <- import_cv_results(
    partial_path = paste0("rf_default._"), #paste0("dnn_","default._"),
    directory = paste0("outputs/",directory)
  )  
  save(cv_models, file = paste0("outputs/",directory,"/","rf_default","_cv.RData"))
  
  lambda_min <- lambda_min(
    pipeline_object = cv_models # the object created from the pipeline object, only works if returnFit=F
  )
  
  final_model_list <- run_pipeline_final(
    feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
    response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
    submission = T,
    method = "glm",
    hyperparam = c("alpha"=1.)#list(c(NULL),c(NULL))#NULL#c("alpha"=0.5)
  )
  
  save(final_model_list, file = paste0("outputs/",directory,"/","glm_","a1.RData"))
  save(lambda_min, file = paste0("outputs/",directory,"/","glm_","a1_lambda.RData"))
  
  # you can "predict(final_model_list[[1]], s = lambda_min[[1]], newx=blablala)" for choosing lambda with optimal cv-score
}
#####################################################






