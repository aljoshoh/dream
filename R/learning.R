### Functions for benchmarking ML algorithms

suppressPackageStartupMessages({
  library(h2o)
  library(dplyr)
  library(tidyverse)
  library(parallel)
  library(glmnet)
  library(glmnetUtils)
  library(rjson)
})

loadRData <- function(
  ### Loads an RData file, and returns it
  ######################################################
  fileName
  ){
  load(fileName)
  get(ls()[ls() != "fileName"])
}######################################################


get_params <- function(
  ### Loads hyperparameters
  ######################################################
  json.file = NULL # json file of paramters, e.g. "params/params.json"
){
  json.file <- as.character(json.file)
  message(paste0('Loading ', json.file))
  params <- fromJSON(json.file)
  
  return(params)
}######################################################


get_features <- function(
  ### Loads features
  ######################################################
  path = NULL # RData file of feature matrix, e.g. "features/features.RData"
){
  path <- as.character(path)
  message(paste0('Loading ', path))
  features <- as.matrix(loadRData(path))
  
  return(features)
}######################################################

dump_features <- function(
  ### Dumps features
  ######################################################
  R.object = NULL, # R obeject of feature matrix
  path = NULL # where to save
){
  message(paste0("Saving ",as.character(deparse(quote(var)))," in ",as.character(path))," ...")
  save(R.object, file = as.character(path))

  return("Done exporting object !")
}######################################################

cv <- function(
  ### Create cross-validation based on row names
  ### Creates a list of length 2 with kfold elements, each 
  ######################################################
  feature_matrix = NULL, # Numerical feature matrix
  phenotype_matrix = NULL, # Numerical pheontype matrix 
  kfold = 5, # How many folds
  seed = 123 # Setting a fixed seed
){
  if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}
  if(!is.matrix(phenotype_matrix)){stop("The phenotype matrix is not a numerical matrix !")}
  message("Check for matched row names...")
  if(!is.null(phenotype_matrix)){
    feature_matrix <- feature_matrix[intersect(row.names(feature_matrix), row.names(phenotype_matrix)),]
    phenotype_matrix <- phenotype_matrix[intersect(row.names(feature_matrix), row.names(phenotype_matrix)),]
  }
  
  message(paste0("Feature matrix is of dimension ",as.character(nrow(feature_matrix)),"x",as.character(ncol(feature_matrix))," !"))
  message(paste0("Phenotype matrix is of dimension ",as.character(nrow(phenotype_matrix)),"x",as.character(ncol(phenotype_matrix))," !"))
  message("Creating feature set names ...")
  
  fold_ids <- cut(1:nrow(feature_matrix), breaks = kfold, labels = F)
  size <- floor(((kfold-1)/kfold)*nrow(feature_matrix))
  message(paste0("Fold size is ",as.character(size)))
  message(paste0("Test size is ", as.character(nrow(feature_matrix)-size)))
  
  set.seed(seed = seed)
  feature_matrix <- feature_matrix[sample(nrow(feature_matrix)),]
  
  test_sets <- lapply(1:kfold, function(x) row.names(feature_matrix[fold_ids == x,]))
  train_sets <- lapply(1:kfold, function(x) row.names(feature_matrix[fold_ids != x,]))
  
  return(list(test_sets = test_sets, 
              train_sets = train_sets
              ))
}######################################################


make_fit <- function(
  ### Runs cross-validation based on row names
  ### Creates a list of length 2 with kfold elements, each 
  ### TODO
  ### -return interesting model parameter
  ######################################################
  feature_matrix = NULL, # Numerical feature matrix
  phenotype_matrix = NULL, # Numerical pheontype matrix 
  folds = NULL, # fold object created in last step
  FUN = function(x){return(x)}, # Pass filtering function for cv sets internally, default is identity function
  hyperparam = c("alpha"=0.5), # Hyperparamter for methods
  method = "glm", # Method
  cvglm = F, # If glm should be cross-validated for the lambda parameter, e.g. for benchmarking
  seed = 123, # Seed for RF method
  returnFit = T # Logical if $param should be returned for all folds at fit, otherwise it returns lambda.min
){
  if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}
  if(!is.matrix(phenotype_matrix)){stop("The feature matrix is not a numerical matrix !")}
  if(length(folds$test_set) != length(folds$train_set)){stop("Folds of testing and training do not match !")}
  
  message("Initialize Method ...")
  if(method %in% c("rf","dnn")){
    h2o.init()
  }
  
  message("Starting to fit ",method," ...")
  score <- matrix(NA,nrow = length(folds$train_set), ncol = ncol(phenotype_matrix))
  param <- as.data.frame(matrix(NA,nrow = length(folds$train_set), ncol = ncol(phenotype_matrix)))
  for(j in 1:ncol(phenotype_matrix)){
    message(paste0("...for drug ", as.character(colnames(phenotype_matrix)[j])))
    for(i in 1:length(folds$train_set)){
      
      message(paste0("...for fold ",as.character(i)))
      y_name <- as.character(colnames(phenotype_matrix)[j])
      y_train <- (phenotype_matrix[folds$train_sets[[i]],j, drop=F])[!is.na(phenotype_matrix[folds$train_sets[[i]],j]),,drop = F]
      x_train <- feature_matrix[names(y_train[,y_name]),]
      
      message(paste0("        response has length ",as.character(length(y_train))," !"))
      x_test <- feature_matrix[folds$test_sets[[i]],]
      y_test <- phenotype_matrix[ folds$test_sets[[i]] ,j]
      
      set.seed(seed)
      #######################METHOD
      if(method == "glm"){
        model <- use_glm(x_train, y_train, x_test, y_test,
                         hyperparam = hyperparam,
                         y_name = y_name,
                         cvglm = cvglm
                         )
      }
      
      if(method == "rf"){
        model <- use_rf(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = y_name,
                        seed = seed
                        )
      }
      
      if(method == "dnn"){
        model <- use_rf(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = as.character(colnames(phenotype_matrix)[j]),
                        seed = seed
        )
      }
      #######################METHOD
      mse <- mean(model$diff*model$diff)
      cor <- cor(model$pred, y_test, use = "complete.obs")
      
      score[i,j] <- cor
      if(returnFit == T){
        param[[i,j]] <- list(model$fit)
      }else{
        param[i,j] <- model$fit$lambda.min
      }
    }
  }
  
  message("End method ...")
  if(method %in% c("rf","dnn")){
    h2o.shutdown(prompt = F)
  }

  return(list(score=score,
              param = param
              ))
}######################################################


use_glm <- function(
  ### Linear Method
  ### TODO
  ### -return interesting model parameter
  ######################################################
  x_train=x_train,
  y_train=y_train,
  x_test=x_test,
  y_test=y_test,
  hyperparam=hyperparam,
  y_name=y_name,
  seed = F,
  cvglm = T
){
  alpha = hyperparam["alpha"]
  y_train <- y_train[,as.character(y_name)]
  
  message("        Fitting...")
  if(cvglm ==T){
    fit <- cv.glmnet(x = x_train, y = y_train, alpha = alpha)
  }else{
    fit <- glmnet::glmnet(x = x_train, y = y_train, alpha = alpha)
  }
  
  message("        Validating...")
  if(cvglm ==T){
    pred <- predict(fit, x_test, s = 'lambda.min')
    diff <- pred - y_test
  }else{
    pred <- NULL
    diff <- NULL
  }
  
  return(list(pred=pred, diff=diff, fit=fit))
}######################################################


use_rf <- function(
  ### Random Forst
  ### TODO
  ### -return interesting model parameter
  ### -tune rf parameter
  ######################################################
  x_train=x_train,
  y_train=y_train,
  x_test=x_test,
  y_test=y_test,
  hyperparam=hyperparam,
  y_name=y_name,
  seed = 123
){
  hyperparams_rf = NULL
  dff <- cbind(x_train,y_train)
  
  message("        Fitting...")
  fit <- h2o.randomForest(x = colnames(x_train), y = y_name, training_frame = as.h2o(dff), seed = seed)
  
  message("        Validating...")
  pred <- predict(fit, as.h2o(x_test))
  pred <- as.data.frame(pred)[,1]
  diff <- pred - y_test
  
  return(list(pred=pred, diff=diff, fit=fit))
}######################################################


use_dnn <- function(
  ### Deep Neural Net
  ### TODO
  ### -return interesting model parameter
  ### -tune rf parameter
  ######################################################
  x_train=x_train,
  y_train=y_train,
  x_test=x_test,
  y_test=y_test,
  hyperparam=hyperparam,
  y_name=y_name,
  seed = 123
){
  hyperparams_rf = NULL
  dff <- cbind(x_train,y_train)
  
  message("        Fitting...")
  fit <- h2o.deeplearning (x = colnames(x_train), y = y_name, training_frame = as.h2o(dff), seed = seed)
  
  message("        Validating...")
  pred <- predict(fit, as.h2o(x_test))
  pred <- as.data.frame(pred)[,1]
  diff <- pred - y_test
  
  return(list(pred=pred, diff=diff, fit=fit))
}######################################################





make_fit_whole <- function(
  ### Makes the model fit on the whole data
  ### Creates a list of length 2 with kfold elements, each 
  ### TODO
  ### -return interesting model parameter
  ######################################################
  feature_matrix = NULL, # Numerical feature matrix
  phenotype_matrix = NULL, # Numerical pheontype matrix 
  FUN = function(x){return(x)}, # Pass filtering function for cv sets internally, default is identity function
  hyperparam = c("alpha"=0.5), # Hyperparamter for methods
  method = "glm", # Method
  seed = 123 # Seed for RF method
){
  if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}
  if(!is.matrix(phenotype_matrix)){stop("The feature matrix is not a numerical matrix !")}
  intersect <- intersect(row.names(feature_matrix), row.names(phenotype_matrix))
  
  message("Initialize Method ...")
  if(method %in% c("rf","dnn")){
    h2o.init()
  }
  
  message("Starting to fit ",method," ...")
  score <- matrix(NA,nrow = 1, ncol = ncol(phenotype_matrix))
  param <- as.data.frame(matrix(NA,nrow = 1, ncol = ncol(phenotype_matrix)))
  for(j in 1:ncol(phenotype_matrix)){
    message(paste0("...for drug ", as.character(colnames(phenotype_matrix)[j])))
    for(i in 1:1){
      
      message(paste0("...for fold ",as.character(i)))
      y_name <- as.character(colnames(phenotype_matrix)[j])
      y_train <- (phenotype_matrix[intersect,j, drop=F])[!is.na(phenotype_matrix[intersect,j]),,drop = F]
      x_train <- feature_matrix[names(y_train[,y_name]),]
      
      message(paste0("        response has length ",as.character(length(y_train))," !"))
      x_test <- NULL #feature_matrix[intersect,]
      y_test <- NULL #phenotype_matrix[ intersect ,j]
      
      set.seed(seed)
      #######################METHOD
      if(method == "glm"){
        model <- use_glm(x_train, y_train, x_test, y_test,
                            hyperparam = hyperparam,
                            y_name = y_name,
                            cvglm = F
        )
      }
      
      if(method == "rf"){
        model <- use_rf(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = y_name,
                        seed = seed
        )
      }
      
      if(method == "dnn"){
        model <- use_rf(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = as.character(colnames(phenotype_matrix)[j]),
                        seed = seed
        )
      }
      #######################METHOD
      mse <- NA
      cor <- NA
      
      score[i,j] <- cor
      param[[i,j]] <- list(model$fit) #<- here we can either asses feature weights, or store other things from model$
    }
  }
  
  message("End method ...")
  if(method %in% c("rf","dnn")){
    h2o.shutdown(prompt = F)
  }
  
  return(list(score=score,
              param = param
  ))
}######################################################
