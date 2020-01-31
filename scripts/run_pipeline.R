### Small scratch on how to use the source codes in R/
### Pipeline on predicting a response matrix based on a feature matrix


run_pipeline_benchmark <- function(
  ### Runs the pipeline and returns the prediction objects
  ### HINT: for more documentation, run View(make_fit)
  ######################################################
  feature_path = NULL, # path to features
  response_path = NULL, # path to response
  submission = T,
  kfold = 10,
  FUN = function(x){return(x)}, # Function on feature set of training data, must return matrix of features (same names as input)
  method = "glm",
  hyperparam = c("alpha"=0.5),
  cvglm = T
){
  if(submission){setwd("storage/groups/cbm01/workspace/dream_aml/")}
  
  # Import objects
  auc <- get_features(feature_path)
  rna <- get_features(response_path)
  
  # Cross-Validation initiation
  CV <- cv(feature_matrix = rna, phenotype_matrix = auc, kfold = kfold, seed = 124)
  
  list_glm <- make_fit(feature_matrix = rna, phenotype_matrix = auc, folds = CV,
                       method = method, 
                       hyperparam = c("alpha"=0.5),
                       cvglm = T, FUN = FUN)
  
  return(list_glm)
}######################################################


### Dump dataframes $rna $auc in the features folder/
  #dump_features(rna, "features/alex_features.RData")
  #dump_features(auc, "features/alex_phenotypes.RData")


# Cross-Validation plot
drug <- 7
plot(log(list_glm$param[1,drug][[1]][[1]]$lambda),list_glm$param[1,drug][[1]][[1]]$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=list_glm$param[4,][[1]][[1]]$name, ylim = c(0,0.01))
points(log(list_glm$param[2,drug][[1]][[1]]$lambda),list_glm$param[2,drug][[1]][[1]]$cvm,pch=19,col="grey")
points(log(list_glm$param[3,drug][[1]][[1]]$lambda),list_glm$param[3,drug][[1]][[1]]$cvm,pch=19,col="blue")
points(log(list_glm$param[4,drug][[1]][[1]]$lambda),list_glm$param[4,drug][[1]][[1]]$cvm,pch=19,col="purple")
points(log(list_glm$param[5,drug][[1]][[1]]$lambda),list_glm$param[5,drug][[1]][[1]]$cvm,pch=19,col="green")
#legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))


## Min lambdas for building final model, make final model with all data, extract weights for the final optimal model
lambda_min <- lapply(1:(dim(list_glm$score)[2]), function(y) median(lapply(1:(dim(list_glm$score)[1]), function(x)list_glm$param[x,y][[1]][[1]]$lambda.min)%>% unlist))

# Cross-Validation plot
drug <- 1
plot(log(list_glm$param[1,drug][[1]][[1]]$lambda),list_glm$param[1,drug][[1]][[1]]$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=list_glm$param[4,][[1]][[1]]$name, ylim = c(0,0.01))
points(log(list_glm$param[2,drug][[1]][[1]]$lambda),list_glm$param[2,drug][[1]][[1]]$cvm,pch=19,col="grey")
points(log(list_glm$param[3,drug][[1]][[1]]$lambda),list_glm$param[3,drug][[1]][[1]]$cvm,pch=19,col="blue")
points(log(list_glm$param[4,drug][[1]][[1]]$lambda),list_glm$param[4,drug][[1]][[1]]$cvm,pch=19,col="purple")
points(log(list_glm$param[5,drug][[1]][[1]]$lambda),list_glm$param[5,drug][[1]][[1]]$cvm,pch=19,col="green")


### Train whole model 
list_glm_whole <- make_fit_whole(feature_matrix = GEXred,
                                 phenotype_matrix = map[,1:5, drop=F],
                                 method = "glm",
                                 hyperparam = c("alpha"=0.3))

fit <- list_glm_whole$param[,drug][[1]][[1]]
beta <- coef(fit, s = lambda_min[[drug]])
beta <- beta[as.logical(beta!=0),, drop=F]
plot(fit, xvar = "lambda"); abline(v=log(lambda_min[[drug]]), col="darkred")






