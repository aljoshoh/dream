#!/usr/bin/env Rscript
setwd("/storage/groups/cbm01/workspace/dream_aml/")
args <- as.numeric(commandArgs(trailingOnly = TRUE)); args <- args[1]
print(paste0("Running with argument: ",as.character(args)))

## Alex' Models
## note that our docker image is using r3.6.1
source("R/learning.R")
source("R/run_pipeline.R")
source("R/select_gene_sc1.R")

models_list <- run_pipeline_benchmark(
  feature_path = "features/alex_features.RData", # path to features
  response_path = "features/alex_phenotypes.RData", # path to response
  submission = T,
  kfold = 10,
  FUN = function(x){return(x)}, # also possible: AnvSigGen 
                                # @phong: the method "make_fit" does not yet return the results of the filtering
  method = "glm",
  hyperparam = c("alpha"=0.5),
  cvglm = T,
  returnFit = F, # if false, then it only returns the lambda
  cvseed = args # supply the parallel processing counter
)

save(models_list, file = paste0("metadata/alex/glm_test_instance",as.character(i),".RData"))





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
