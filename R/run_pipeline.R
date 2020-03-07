### Pipeline functions on predicting a response matrix based on a feature matrix


run_pipeline_benchmark <- function(
  ### Runs the pipeline and returns the prediction objects
  ### HINT: for more documentation, run View(make_fit)
  ######################################################
  feature_path = NULL, # path to features
  response_path = NULL, # path to response
  submission = T,
  kfold = 10,
  FUN = function(x,y){res<-lapply(1:ncol(y),function(z) colnames(x));return(res)}, # Function on feature set of training data, must return matrix of features (same names as input)
  method = "glm",
  hyperparam = c("alpha"=0.5),
  cvglm = T,
  returnFit = F, # if false, then it only returns the lambda
  cvseed = 123,
  CVBuilt = NULL,
  stack = F,
  args = NULL# <---------------------------------------------------------------------------------------------------------------- needs to be implemented consistently
){
  if(submission){setwd("/storage/groups/cbm01/workspace/dream_aml/")}
  
  # Import objects
  rna <- get_features(feature_path)
  auc <- get_features(response_path)
  
  if(stack){
    rna <- get_features(feature_path, matrixfy = F)
  }
  
  # Cross-Validation initiation
  if(!is.null(CVBuilt)){
    message("Using CV-Built !")
    CV <- CVBuilt
  } else {
    CV <- cv(feature_matrix = rna, phenotype_matrix = auc, kfold = kfold, seed = cvseed)
  }
  
  list_glm <- make_fit(feature_matrix = rna, phenotype_matrix = auc, folds = CV,
                       method = method, 
                       hyperparam = hyperparam,
                       cvglm = T, FUN = FUN, returnFit = returnFit, stack = stack, args = args)
  
  return(list_glm)
}######################################################


### Dump dataframes $rna $auc in the features folder/
  #dump_features(rna, "features/alex_features.RData")
  #dump_features(auc, "features/alex_phenotypes.RData")


plot_cv <- function(
  ### Plot the cross-validations for each fold
  ### Only works if returnFit=F
  ######################################################
  drug_index = 1, # which drug index to plot
  pipeline_object = NULL # the object created from the pipeline object
){
  if(pipeline_object$returnFit){stop("This function only works if only lambda is returned !!!")}
  drug <- drug_index
  colors <- rainbow(dim(pipeline_object$score[1]))
  
  plot(log(list_glm$param[1,drug][[1]][[1]]$lambda),list_glm$param[1,drug][[1]][[1]]$cvm,pch=19,col=colors[1],xlab="log(Lambda)",ylab=list_glm$param[4,][[1]][[1]]$name)
  for(i in 2:(dim(object$score[1]))){
    points(log(list_glm$param[i,drug][[1]][[1]]$lambda),list_glm$param[i,drug][[1]][[1]]$cvm,pch=19,col=colors[i])
  }
  
  return(paste0("Plot of benchmark for index: ",as.character(drug_index),"..."))
}######################################################



lambda_min <- function(
  ### get minimum lambdas for each drug in the pipeline object
  ### Only works if returnFit=F
  ### TODO
  ### Test this function
  ######################################################
  pipeline_object = NULL # the object created from the pipeline object
){
  #if(pipeline_object$returnFit){stop("This function only works if only lambda is returned !!!")}
  lambda_min <- lapply(1:(dim(pipeline_object$score)[2]), 
                       function(y) median(lapply(1:(dim(pipeline_object$score)[1]), 
                                                 function(x) pipeline_object$param[x,y][[1]][[1]]$lambda.min)%>% unlist))
  return(lambda_min)
}######################################################



run_pipeline_final <- function(
  ### Builds the final model based on all the data
  ######################################################
  feature_path = NULL, # path to features
  response_path = NULL, # path to response
  submission = T,
  FUN = function(x){return(x)}, # Function on feature set of training data, must return matrix of features (same names as input)
  method = "glm",
  hyperparam = c("alpha"=0.5),
  stack = F
){
  if(submission){setwd("/storage/groups/cbm01/workspace/dream_aml/")}
  
  # Import objects
  rna <- get_features(feature_path)
  auc <- get_features(response_path)
  
  ### Train whole model
  list_glm_whole <- make_fit_whole(feature_matrix = rna,
                                   phenotype_matrix = auc,
                                   method = method,
                                   hyperparam = hyperparam)

  fit <- list()
  for(drug in 1:ncol(auc)){
    fit[[drug]] <- list_glm_whole$param[,drug][[1]][[1]]
  }
    #Potentially return the betas here
    #beta <- coef(fit[[drug]], s = lambda_min[[drug]])
    #beta <- beta[as.logical(beta!=0),, drop=F]
    #plot(fit, xvar = "lambda"); abline(v=log(lambda_min[[drug]]), col="darkred")
  
  
  return(list_glm_whole)
}######################################################




##############
## DEPRECATED
##############


# Cross-Validation plot
#drug <- 7

#points(log(list_glm$param[2,drug][[1]][[1]]$lambda),list_glm$param[2,drug][[1]][[1]]$cvm,pch=19,col="grey")
#points(log(list_glm$param[3,drug][[1]][[1]]$lambda),list_glm$param[3,drug][[1]][[1]]$cvm,pch=19,col="blue")
#points(log(list_glm$param[4,drug][[1]][[1]]$lambda),list_glm$param[4,drug][[1]][[1]]$cvm,pch=19,col="purple")
#points(log(list_glm$param[5,drug][[1]][[1]]$lambda),list_glm$param[5,drug][[1]][[1]]$cvm,pch=19,col="green")
#legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))


## Min lambdas for building final model, make final model with all data, extract weights for the final optimal model


# Cross-Validation plot
#drug <- 1
#plot(log(list_glm$param[1,drug][[1]][[1]]$lambda),list_glm$param[1,drug][[1]][[1]]$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=list_glm$param[4,][[1]][[1]]$name, ylim = c(0,0.01))
#points(log(list_glm$param[2,drug][[1]][[1]]$lambda),list_glm$param[2,drug][[1]][[1]]$cvm,pch=19,col="grey")
#points(log(list_glm$param[3,drug][[1]][[1]]$lambda),list_glm$param[3,drug][[1]][[1]]$cvm,pch=19,col="blue")
#points(log(list_glm$param[4,drug][[1]][[1]]$lambda),list_glm$param[4,drug][[1]][[1]]$cvm,pch=19,col="purple")
#points(log(list_glm$param[5,drug][[1]][[1]]$lambda),list_glm$param[5,drug][[1]][[1]]$cvm,pch=19,col="green")


