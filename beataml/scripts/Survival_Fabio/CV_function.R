RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

crossvalidation <- function(entity, nFold=5, random=T, seed=F) {
  if(!random){
    if(!seed){
      set.seed(754)
    } else {
      set.seed(seed)
    }
  }
  
  vec <- entity[sample(length(entity))]
  
  size <- floor(length(vec)/nFold)
  binSize <- rep(size, nFold)
  modulo <- length(vec) %% size
  if (modulo > 0) {
    binSize[1:modulo] <- binSize[1:modulo] +1
  }
  
  binEndIdx <- cumsum(binSize)
  binStardIdx <- c(1, binEndIdx[1:(nFold-1)]+1)
  
  cv <- list()
  for (i in 1:nFold){
    
    xTrainIdx <- binStardIdx[(i %% nFold) + 1]:binEndIdx[(i %% nFold) + 1]
    xTrain <- vec[xTrainIdx]
    
    train <- vec[setdiff(1:length(vec), xTrainIdx)]
    
    cv[[paste("fold_", i, sep="")]] <- list(xTrain=xTrain, train=train)
  }
  return(cv) 
}

RunCV <- function(Combination.Hyperparameters, InputExpression, Survival) {
  for (i in 1:nrow(Combination.Hyperparameters) ){
    
    cv <- crossvalidation(colnames(Survival), 5)
    
    Rp <- list()
    lambda <- list()
    pred <- c()
    for (j in 1:length(cv)) {# for each cv set
      temp.xTrain <- InputExpression[, cv[[j]]$train]
      temp.xPred <- InputExpression[,cv[[j]]$xTrain]
      temp.yTrain <- Survival[,cv[[j]]$train]
      temp.yPred <-Survival[,cv[[j]]$xTrain,]
      
      tempTrain <- as.data.frame(t(rbind(temp.xTrain, temp.yTrain)))
      tempPred <-as.data.frame(t(rbind(temp.xPred, temp.yPred)))
      # 1) train model
      fit <- rfsrc(Surv(overallSurvival, vitalStatus) ~ .,data = tempTrain, ntree=Combination.Hyperparameters[i,1], mtry =Combination.Hyperparameters[i,2], splitrule ='logrank')
      # 2) predict on validation set
      survival.results <- predict.rfsrc(fit, newdata = tempPred)
      # 3) calculate performance on validaiton set
      Rp[[j]] <- survival.results$err.rate
      #Rp[[i]] <- apply(xTrain_pred, 2, function(x) if (sd(x)>0) {score_metric(x, drug[cv[[i]]$xTrain])} else{NA})
    }
    
    nMax <- max(unlist(lapply(Rp, length))) # maximal pathlength occuring in all the cv regulatization paths
    Rp <- lapply(Rp, function(x) c(x, rep(NA, nMax-length(x))))
    Rp <- matrix(unlist(Rp), nrow=length(Rp), byrow=T, dimnames = list(NULL,names(Rp[[1]])))# fill up NAs and make a full matrix
    
    # 5) chose representative
    Combination.Hyperparameters[i,"ErrRate" ] <- median(Rp, na.rm = TRUE)
  }
  BestCombination <- Combination.Hyperparameters[which.min(Combination.Hyperparameters$ErrRate),]
  return(BestCombination)
}
