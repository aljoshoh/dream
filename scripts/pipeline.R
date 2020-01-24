dump_features(rna, "features/alex_features.RData")
dump_features(auc, "features/alex_phenotypes.RData")

auc <- get_features("features/alex_phenotypes.RData")
rna <- get_features("features/alex_features.RData")


CV <- cv(feature_matrix = rna, phenotype_matrix = auc, kfold = 5, seed = 123)

list_glm <- make_fit(feature_matrix = rna[,1:12], phenotype_matrix = auc[,1:1,drop=F], folds = CV,
                 method = "dnn", 
                 hyperparam = c("alpha"=0.5))



# Cross-Validation plot
plot(log(list_glm$param[1,][[1]][[1]]$lambda),list_glm$param[1,][[1]][[1]]$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=list_glm$param[4,][[1]][[1]]$name, ylim = c(1000,2200))
points(log(list_glm$param[2,][[1]][[1]]$lambda),list_glm$param[2,][[1]][[1]]$cvm,pch=19,col="grey")
points(log(list_glm$param[3,][[1]][[1]]$lambda),list_glm$param[3,][[1]][[1]]$cvm,pch=19,col="blue")
points(log(list_glm$param[4,][[1]][[1]]$lambda),list_glm$param[4,][[1]][[1]]$cvm,pch=19,col="purple")
points(log(list_glm$param[5,][[1]][[1]]$lambda),list_glm$param[5,][[1]][[1]]$cvm,pch=19,col="green")
#legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))
plot(list_glm$param[2,][[1]][[1]])






h2o.init()
df <- cbind(rna[1:171,1:10],auc[1:171,1,drop=F])
test <- h2o.randomForest(x = colnames(rna[1:171,1:10]), y = colnames(auc[1:171,1,drop=F]), training_frame = as.h2o(df),
                         seed = 123)

pred <- predict(test, as.h2o(rna[172:213,1:10]))
cor(as.data.frame(pred)[,1], auc[172:213,1], use="complete.obs")^2
#rss <- sum((as.data.frame(pred)[,1] - auc[172:213,1]) ^ 2, na.rm=T)  ## residual sum of squares
#tss <- sum((auc[172:213,1] - mean(auc[172:213,1], na.rm=T)) ^ 2, na.rm=T)  ## total sum of squares
#rsq <- 1 - rss/tss; rsq
#mean((as.data.frame(pred)[,1]- auc[172:213,1])^2, na.rm= T)
#h2o.performance(model = test, newdata = as.h2o(cbind(rna[172:213,1:10],auc[172:213,1,drop=F])), valid = F)
h2o.shutdown()



