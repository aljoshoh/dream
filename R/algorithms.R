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