library(randomForestSRC)
library(survival)

ExpressionMatrix <- data.frame(read.table('~/rnaseq.csv', header = T, sep = ','))
rownames(ExpressionMatrix) <- ExpressionMatrix[,1]
ExpressionMatrix <- ExpressionMatrix[,-c(1,2)]
colnames(ExpressionMatrix) <- gsub(".", "-",colnames(ExpressionMatrix), fixed =T)
colnames(ExpressionMatrix) <- gsub("X", "",colnames(ExpressionMatrix), fixed =T)

Survival <- data.frame(read.table('~/response.csv', header = T, sep = ',', row.names = 1))
ExpressionMatrix <- ExpressionMatrix[, colnames(ExpressionMatrix )%in% rownames(Survival)]
Survival <- Survival[ colnames(ExpressionMatrix ),]
Survival$vitalStatus <- as.numeric(factor(Survival$vitalStatus))
for (i in 1:dim(Survival)[1]) {
  if(Survival$vitalStatus[i] == "2"){
    Survival$vitalStatus[i]<- 1
    } else {
      Survival$vitalStatus[i]<- 0
  }
}
Survival <- as.data.frame(t(Survival))
colnames(Survival) <- colnames(ExpressionMatrix)

# GeneSet.Metzeleret <- as.data.frame(read.csv("~/Desktop/DREAMChallenge_BeatAML/Survival/SurvivalAssociatedGenes_MetzeleretAl.csv", sep=";"))
# GeneSet.Metzeleret <- as.data.frame(GeneSet.Metzeleret[,1])
# GeneSet.Nagyet1 <-as.data.frame(read.csv("~/Desktop/DREAMChallenge_BeatAML/Survival/SurvivalAssociatedgenes_NagyetAl.csv", sep=";"))
# GeneSet.Nagyet1 <- as.data.frame(GeneSet.Nagyet1[!(GeneSet.Nagyet1[,1]%in% ""),  1])
# GeneSet.Nagyet2 <-as.data.frame(read.csv("~/Desktop/DREAMChallenge_BeatAML/Survival/SurvivalAssociatedgenes_NagyetAl_2.csv", sep=";"))
# GeneSet.Nagyet2 <- as.data.frame(GeneSet.Nagyet2[!(GeneSet.Nagyet2[,1]%in% ""),  1])

ExpressionMatrix <- ExpressionMatrix[rowSums(ExpressionMatrix) >5,]
# ExpressionMatrix.Metzeleret <- ExpressionMatrix[GeneSet.Metzeleret[,1], ]
# ExpressionMatrix.Nagyet1 <- ExpressionMatrix[GeneSet.Nagyet1[,1], ]
# ExpressionMatrix.Nagyet2 <- ExpressionMatrix[GeneSet.Nagyet2[,1], ]

###### Random forest for survival

source("~/Survival_Fabio/CV_function.R")

ntree <- seq(from=10, to= 5000, by = 3000)
mtry <- seq(from= 10, to= 2000, by= 3000)

Combination.Hyperparameters <- as.data.frame(as.matrix(tidyr::crossing(ntree, mtry)))
InputExpression <- ExpressionMatrix

Best.Param <- RunCV(Combination.Hyperparameters, InputExpression, Survival)

# Retrain model on entire dataset usign optimal parameters
InputGeneral <- as.data.frame(t(rbind(ExpressionMatrix, Survival)))
fit.Best <- rfsrc(Surv(overallSurvival, vitalStatus) ~ .,data = InputGeneral, ntree=Best.Param[1,1], mtry =Best.Param[1,2], splitrule ='logrank')

