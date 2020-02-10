## Random forest model stacking

library(randomForest)
library(caret)
source("R/learning.R")
source("R/run_pipeline.R")
source("R/select_gene_sc1.R")
source("R/algorithms.R")


feature_path = "features/alex_features.RData" # path to features
response_path = "features/alex_phenotypes.RData" # path to response
rna <- get_features(feature_path)
auc <- get_features(response_path)

rna <- rna[,1:1000]
auc <- auc[,1:2]

dump_features(auc, path = "features/alex_phenotypes_red.RData")
dump_features(rna, path = "features/alex_features_red.RData")

test <- randomForest(x = rna[!is.na(auc[,1]),], 
                     y = (auc[,1])[!is.na(auc[,1])],
                     
                     )

plot(predict(test, newdata = rna[!is.na(auc[,1]),]), (auc[,1])[!is.na(auc[,1])])




seed <- 7
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



# train model
control <- trainControl(method="repeatedcv", number=10, repeats=1)
tunegrid <- expand.grid(.mtry=c(15,200), .ntree= 2000) # c(1000, 1500, 2000, 2500) 1:15

custom <- train(x = rna[!is.na(auc[,1]),],
                y = (auc[,1])[!is.na(auc[,1])],
                method=customRF, 
                metric="RMSE", 
                tuneGrid=tunegrid, 
                trControl=control)
summary(custom)
plot(custom)




plot(custom)