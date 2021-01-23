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
library(survcomp)
library(tidyr)
source("R/algorithms.R")
source("R/general.R")

#### how to run on cluster
# ./SUBMIT.sh PIPELINE.R numberofargs 1
args <- args[1] #<- new command 
print(paste0("Running with argument: ",as.character(args)))
##########################


mtry <- seq(from= 1, to= 8, by= 1) #rna: c(10,20,40,70,seq(from= 100, to= 1000, by= 100))
ntree <- c(100, 200, 500, 1000)
grid <- expand_grid(mtry, ntree)%>%as.data.frame
#param <- list(mtry,ntree) 
param <- lapply(1:nrow(grid), function(x) list(grid[x,1],grid[x,2]))[[args]]

here<- "cox"
modelsa <- get_features(paste0("outputs/mut-surv/",here,"_default_cv.RData"), matrixfy = F)
modelsb <- get_features(paste0("outputs/rna-surv/",here,"_default_cv.RData"), matrixfy = F)
modelsc <- get_features(paste0("outputs/clin-surv/",here,"_default_cv.RData"), matrixfy = F)
modelsd <- get_features(paste0("outputs/auc-surv/",here,"_default_cv.RData"), matrixfy = F)
here<- "rfsurv"
modelse <- get_features(paste0("outputs/mut-surv/",here,"_default_cv.RData"), matrixfy = F)
modelsf <- get_features(paste0("outputs/rna-surv/",here,"_default_cv.RData"), matrixfy = F)
modelsg <- get_features(paste0("outputs/clin-surv/",here,"_default_cv.RData"), matrixfy = F)
modelsh <- get_features(paste0("outputs/auc-surv/",here,"_default_cv.RData"), matrixfy = F)

#h2o.init()
#modelsc$param_paths <- modelsc$param
#modelsf$param_paths <- modelsc$param
#modelsi$param_paths <- modelsc$param

### SCRIPT PARAMETER
directory <- "auc-surv"#"mut" #"rna"
descriptor <- here # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response_1.RData") # path to response
auc <- get_features(feature_path)
surv_auc <- get_features(response_path)

### SCRIPT PARAMETER
directory <- "rna-surv"#"mut" #"rna"
descriptor <- here # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response_1.RData") # path to response
rna <- get_features(feature_path)
surv_rna <- get_features(response_path)

### SCRIPT PARAMETER
directory <- "mut-surv"#"mut" #"rna"
descriptor <- here # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response_1.RData") # path to response
mut <- get_features(feature_path)
surv_mut <- get_features(response_path)

### SCRIPT PARAMETER
directory <- "clin-surv"#"mut" #"rna"
descriptor <- here # the descriptor does not care here
####################
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response_1.RData") # path to response
clin <- get_features(feature_path)
surv_clin <- get_features(response_path)

if(FALSE){
    preds_stack <- list()
    for(y in 1:1){ 
      
      model1 <- unlist(lapply(1:length(modelsa$cv$test_sets), function(x){ tmp<-predict(modelsa$param[[x,y]][[1]], newx = as.matrix(mut[modelsa$cv$test_sets[[x]],modelsa$gene_names_filtered[[x,y]][[1]]]),s = 'lambda.min' )[,1]; names(tmp)=row.names(mut[modelsa$cv$test_sets[[x]],]); return(tmp)}))
      model2 <- unlist(lapply(1:length(modelsb$cv$test_sets), function(x){ tmp<-predict(modelsb$param[[x,y]][[1]], newx = as.matrix(rna[modelsb$cv$test_sets[[x]],modelsb$gene_names_filtered[[x,y]][[1]]]),s = 'lambda.min')[,1]; names(tmp)=row.names(rna[modelsb$cv$test_sets[[x]],]); return(tmp)}))
      model3 <- unlist(lapply(1:length(modelsc$cv$test_sets), function(x){ tmp<-predict(modelsc$param[[x,y]][[1]], newx = as.matrix(clin[modelsc$cv$test_sets[[x]],modelsc$gene_names_filtered[[x,y]][[1]]]),s = 'lambda.min')[,1]; names(tmp)=row.names(clin[modelsc$cv$test_sets[[x]],]); return(tmp)}))
      model4 <- unlist(lapply(1:length(modelsd$cv$test_sets), function(x){ tmp<-predict(modelsd$param[[x,y]][[1]], newx = as.matrix(auc[modelsd$cv$test_sets[[x]],modelsd$gene_names_filtered[[x,y]][[1]]]),s = 'lambda.min')[,1]; names(tmp)=row.names(auc[modelsd$cv$test_sets[[x]],]); return(tmp)}))
      
      model5 <- unlist(lapply(1:length(modelse$cv$test_sets), function(x){ tmp<-predict.rfsrc(modelse$param[[x,y]][[1]], newdata = as.data.frame(mut[modelse$cv$test_sets[[x]],]))$predicted; names(tmp)=row.names(mut[modelse$cv$test_sets[[x]],]); return(tmp)}))
      model6 <- unlist(lapply(1:length(modelsf$cv$test_sets), function(x){ tmp<-predict.rfsrc(modelsf$param[[x,y]][[1]], newdata = as.data.frame(rna[modelsf$cv$test_sets[[x]],]))$predicted; names(tmp)=row.names(rna[modelsf$cv$test_sets[[x]],]); return(tmp)}))
      model7 <- unlist(lapply(1:length(modelsg$cv$test_sets), function(x){ tmp<-predict.rfsrc(modelsg$param[[x,y]][[1]], newdata = as.data.frame(clin[modelsg$cv$test_sets[[x]],]))$predicted; names(tmp)=row.names(clin[modelsg$cv$test_sets[[x]],]); return(tmp)}))
      model8 <- unlist(lapply(1:length(modelsh$cv$test_sets), function(x){ tmp<-predict.rfsrc(modelsh$param[[x,y]][[1]], newdata = as.data.frame(auc[modelsh$cv$test_sets[[x]],]))$predicted; names(tmp)=row.names(auc[modelsh$cv$test_sets[[x]],]); return(tmp)}))
      intsec <- intersect(intersect(intersect(names(model1),names(model2)), names(model3)), names(model4))
      temp <- cbind(model1[intsec], 
                    model2[intsec], 
                    model3[intsec],
                    model4[intsec],
                    model5[intsec],
                    model6[intsec],
                    model7[intsec],
                    model8[intsec] 
                    );
      preds_stack[[y]] <- temp
      print(paste0(as.character(y/ncol(surv_rna)*100),"%"))
    }
    
    stack_features <- preds_stack
    for(i in 1:length(stack_features)){
      colnames(stack_features[[i]]) = c("mut.cox","rna.cox","clin.cox","auc.cox","mut.rf","rna.rf","clin.rf","auc.rf")
    }
    dump_features(stack_features, path = "features/stacked_models_sc2_features_update.RData")
    surv_stack <- surv_rna[intsec,]
    dump_features(surv_stack, path = "features/stacked_models_sc2_response_update.RData")

}

models_list_stacked <- run_pipeline_benchmark(
  feature_path = "features/stacked_models_sc2_features_update.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
  response_path = "features/stacked_models_sc2_response_update.RData", # path to response
  submission = T,
  kfold = 10, 
  method = "rfsurv",
  hyperparam = param, #list(c(333),c(500)), # c("alpha"=0.5),
  cvglm = T,
  returnFit = T, # if false, then it only returns the lambda
  cvseed = NULL, # using the supplied CV built
  #CVBuilt = modelsa$cv,
  stack = T
)
models_list_stacked$score[,1] %>% replace_na(0.5) %>% mean


save(models_list_stacked, file = paste0("outputs/stacked_models_sc2","_cv_update_",as.character(args),"_",as.character(round(mean(models_list_stacked$score[,1]%>%replace_na(0.5)),3)),".RData"))

if(FALSE){
models_list_stacked <- run_pipeline_final(
  feature_path = "features/stacked_models_sc2_features_update.RData", # path to features, this time as list orderer like the drugs in the response path file !!!
  response_path = "features/stacked_models_sc2_response_update.RData", # path to response
  submission = T,
  method = "rfsurv",
  hyperparam = list(c(1),c(500)), #c("alpha"=0.5), #list(c(333),c(500)), # c("alpha"=0.5),
  stack = T
)
}
save(models_list_stacked, file = "outputs/stacked_models_sc2.RData")
