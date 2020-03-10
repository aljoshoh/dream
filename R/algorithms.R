### ML algorithms for classical x_train,y_train,x_test,y_test strucuture


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


use_cox <- function(
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
  x_train <- data.matrix(as.data.frame(x_train))
  message("        Fitting...")
  if(cvglm ==T){
    fit <- cv.glmnet(x = x_train, y = y_train, alpha = alpha, family="cox")
  }else{
    fit <- glmnet::glmnet(x = x_train, y = y_train, alpha = alpha, family = "cox")
  }
  
  message("        Validating...")
  if(cvglm ==T){
    x_train <<- x_train
    x_test <- data.matrix(as.data.frame(x_test))
    pred <- predict(fit, x_test, s = 'lambda.min')
    diff <- NULL # does not work for survival
  }else{
    pred <- NULL
    diff <- NULL
  }
  
  return(list(pred=pred, diff=diff, fit=fit))
}######################################################

use_rfSurvival <- function(x_train=x_train,
y_train=y_train,
x_test=x_test,
y_test=y_test,
hyperparam=hyperparam,
y_name=y_name,
seed = F,
) {
  Train <- as.data.frame(t(rbind(x_train, y_train)))
  Pred <-as.data.frame(t(rbind(x_test, y_test)))
  
  # 1) train model
  fit <- rfsrc(Surv(overallSurvival, vitalStatus) ~ .,data = Train, ntree=hyperparam$mtree, mtry =hyperparam$mtry, splitrule ='logrank')
  # 2) predict on validation set
  survival.results <- predict.rfsrc(fit, newdata = Pred)
  # 3) calculate performance on validaiton set
  Error <- survival.results$err.rate
  
  return (list(pred = survival.results, diff = NULL, fit = fit))
}




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
  seed = NULL,
  kfold = NULL
){
  # Initialize method
  customRF <- list(type = "Regression", library = "randomForest", loop = NULL)
  customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
  customRF$grid <- function(x, y, len = NULL, search = "grid") {}
  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
  }
  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata = newdata
    )
  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")
  customRF$sort <- function(x) x[order(x[,1]),]
  customRF$levels <- function(x) x$classes
  set.seed(seed)
  
  message("        Fitting...")
  
  if(is.null(kfold)){ # if whole data training
    control <- trainControl(method="none")
  }else{
    control <- trainControl(method="repeatedcv", number=kfold, repeats=1)
  }
  
  if(is.null(hyperparam[[1]]) & is.null(hyperparam[[2]])){
    tunegrid <- expand.grid(.mtry=round(ncol(x_train)/3), .ntree= 500) # here the hyperparams are feeded
  }else{
    tunegrid <- expand.grid(.mtry=hyperparam[[1]], .ntree= hyperparam[[2]]) # here the hyperparams are feeded
  }
  
  fit <- train(x = as.data.frame(x_train),
                  y = as.numeric(y_train),
                  method=customRF, 
                  metric="RMSE", 
                  tuneGrid=tunegrid, 
                  trControl=control)
  
  
  message("        Validating...")
  pred <- predict(fit, x_test)
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
  fit <<- h2o.deeplearning (x = colnames(x_train), y = y_name, training_frame = as.h2o(dff), seed = seed)
  
  message("        Validating...")
  if(!is.null(x_test)){
    x_test <- as.h2o(x_test)
    pred <- predict(fit, x_test)
    pred <- as.data.frame(pred)[,1]
  }else{
    pred <- NULL
  }
  
  diff <- pred - y_test
  fit <- h2o.saveModel(object=fit, path="metadata/h2odnn/", force=TRUE)
  
  return(list(pred=pred, diff=diff, fit=fit))
}######################################################
