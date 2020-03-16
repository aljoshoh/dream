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
# ./SUBMIT.sh PIPELINE.R numberofargs 1
#args <- args[1] #<- new command 
print(paste0("Running with argument: ",as.character(args)))
##########################

### SCRIPT PARAMETER
directory <- "rna-surv"#"mut" #"rna"
descriptor <- "rfsurv" # the descriptor means the method in this script, not the same as in PREPROCESS.R

mtry <- seq(from= 1, to= 1000, by= 100) #rna:seq(from= 1, to= 1000, by= 50) 
ntree <- seq(from=100, to= 1000, by = 100) 
combos <- expand.grid(mtry, ntree) %>% t
combos_local <- cut_df(df = combos, totaln = 40, index = args)
print(ncol(combos_local))

#clin:param=list(c(10),c(500))
#mut:param=list(c(16),c(100))
#auc:param=list(c(),)
####################
if(descriptor=="dnn"){h2o.init(port=8508)}
cvscore <- c()
for(i in 1:nrow(combos_local)){
  param <- list(c(combos_local[i,1]),c(combos_local[i,2])) 
  models_list <- run_pipeline_benchmark(
    feature_path = paste0("features/",directory,"/",descriptor,"_features.RData"), # path to features
    response_path = paste0("features/",directory,"/",descriptor,"_response_",as.character(1),".RData"), # path to response, all the same for survival
    submission = T,
    kfold = 10, 
    method = c(descriptor),
    hyperparam = param,
    cvglm = T,
    returnFit = T, # if false, then it only returns the lambda
    cvseed = 1, #CV-seed for tuning was set to 1, setting 2 for creating the stacking features !
    FUN = AnvSigSurvFeature,
    args = args
  )
  cvscore[i] <- mean(unlist(models_list$score), na.rm=T)
  print(i/nrow(combos)*100)
}
combos_local$cvscore <- cvscore
save(combos_local, file = paste0("outputs/",directory,"/tuning/",descriptor,"_performance_10fold_cvseed1_instance",as.character(args),"100var.RData"))
#save(models_list, file = paste0("outputs/",directory,"/",descriptor,"_default._10fold_cvseed1_instance",as.character(args),".RData"))

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
    hyperparam = list(c(NULL),c(NULL)),
    FUN = AnvSigSurvFeature #AnvSigSurvFeature
  )
  
  save(final_model_list, file = paste0("outputs/",directory,"/",descriptor,"_default.RData"))
  
  
  # you can "predict(final_model_list[[1]], s = lambda_min[[1]], newx=blablala)" for choosing lambda with optimal cv-score
}
#####################################################





########## Gridsearch
if(FALSE){
  library(reshape2)
  s <- spread(combos_local, key = "Var1", value = cvscore)
  row.names(s) = s$Var2
  s <- s[,-1]
  s <- as.matrix(s)
  s <- melt(s)
  ggplot(s, aes(x = Var2, y = Var1)) + 
         geom_raster(aes(fill=value)) + 
         scale_fill_gradient(low="grey90", high="red") +
         labs(x="letters", y="LETTERS", title="Matrix") +
         theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                                                     axis.text.y=element_text(size=9),
                                                     plot.title=element_text(size=11))
}



