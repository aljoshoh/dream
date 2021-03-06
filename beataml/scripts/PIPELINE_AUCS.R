#!/usr/bin/env Rscript
setwd("/storage/groups/cbm01/workspace/dream_aml/")
args <- as.numeric(commandArgs(trailingOnly = TRUE)); args <- args # get the submission script argument
print("Arguments:")
print(args)
library(survcomp)
#args <- 1 # TEMPORALY LOCALLY NEEDED FOR TESTING

## Alex' Models
## note that our docker image is using r3.6.1
source("R/learning.R")
source("R/run_pipeline.R")
source("R/feature_selection.R")
source("R/algorithms.R")
source("R/general.R")

#### how to run on cluster
# ./SUBMIT.sh PIPELINE_AUCS.R numberofargs 1
args <- args[1] #<- new command 
print(paste0("Running with argument: ",as.character(args)))
##########################

### SCRIPT PARAMETER
directory <- "rna-auc"#"mut" #"rna"
descriptor <- "rf" # the descriptor means the method in this script, not the same as in PREPROCESS.R
param <- list(c(NULL),c(NULL)) #list(c(333),c(500)) # c("alpha"=1.)

mtry <- c(seq(from= 1, to= 10, by= 3),20,27,35,40) #rna: c(10,20,40,70,seq(from= 100, to= 1000, by= 100))
ntree <- c(100, 200, 500, 1000)
param <- list(mtry,ntree) 

####################
if(descriptor=="dnn"){h2o.init(port=8508)}
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
  FUN = AnvSigNumFeature,
  args = args
)
models_list$cvscore <- mean(unlist(models_list$score), na.rm=T)
save(models_list, file = paste0("outputs/",directory,"/",descriptor,"_default._10fold_cvseed1_instance",as.character(args),".RData"))


h2o.shutdown()


########### FURTHER STEPS FOR BUILDING THE FULL MODEL
# You can do these things also in a separate scripts, since these are not that computationally heavy
if(FALSE){

  cv_models <- import_cv_results(
    partial_path = paste0(descriptor,"_default._"),
    directory = paste0("outputs/",directory)
  )  
  #cv_models <- import_cv_results(
  #  +     partial_path = paste0(descriptor,"_cv_update_"),
  #  +     directory = paste0("outputs/",directory)
  #  + )  
  save(cv_models, file = paste0("outputs/",directory,"/",descriptor,"_default_cv.RData"))
  
  if(descriptor=="glm"){
    lambda.min <- lambda_min(
      pipeline_object = cv_models # the object created from the pipeline object, only works if returnFit=F
    )
    save(lambda.min, file =  paste0("outputs/",directory,"/",descriptor,"_default_lambda.min.RData"))
  }
  
  model_object  <- loadRData(paste0("outputs/",directory,"/",descriptor,"_default_cv.RData"))
  if(descriptor=="rf"){
    hypermin <- lapply(1:122,
          function(x)
            list(
            median(unlist(lapply(1:10, function(y) model_object$param[[y,x]][[1]]$bestTune$mtry))),
            median(unlist(lapply(1:10, function(y) model_object$param[[y,x]][[1]]$bestTune$ntree)))
            )
          )
    save(hypermin, file =  paste0("outputs/",directory,"/",descriptor,"_default_hyper.min.RData"))
  }
  
  ### Train model
  if(descriptor=="dnn"){h2o.init(port=8510)}
  final_model_list <- run_pipeline_final(
    feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
    response_path = paste0("features/",directory,"/",descriptor,"_response.RData"), # path to response
    submission = T,
    method = descriptor,
    hyperparam = hypermin,
    FUN = AnvSigNumFeature
  )
  
  save(final_model_list, file = paste0("outputs/",directory,"/",descriptor,"_default.RData"))
  
  
  # you can "predict(final_model_list[[1]], s = lambda_min[[1]], newx=blablala)" for choosing lambda with optimal cv-score
}
#####################################################






