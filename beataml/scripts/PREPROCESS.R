numberofargs <- 1 # if sequential, set to 1
###############

directory <- "auc-surv"#"mut" #"rna"
descriptor <- directory
feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response

rna <- get_features(feature_path)
auc <- get_features(response_path)
rna <- rna[intersect(row.names(rna),row.names(auc)),]
auc <- auc[intersect(row.names(rna),row.names(auc)),]
print(dim(rna))
print(dim(auc))
descriptor <- "cox" #method #"coxrf"
dump_features(rna, path = paste0("features/",directory,"/",descriptor,"_features.RData"))
dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response.RData"))

auc_pre <- auc
for(args in 1:numberofargs){
  if(numberofargs ==1){
    auc <- auc_pre
  } else {
    auc <- cut_df(auc_pre, numberofargs,args)
  }
  dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response_",as.character(args),".RData"))
}




if(FALSE){
  numberofargs <- 40 # if sequential, set to 1
  ###############
  
  directory <- "clin-auc"#"mut" #"rna"
  descriptor <- directory
  feature_path = paste0("features/",directory,"/",descriptor,"_features.RData") # path to features
  response_path = paste0("features/",directory,"/",descriptor,"_response.RData") # path to response
  
  rna <- get_features(feature_path)
  auc <- get_features(response_path)
  rna <- rna[intersect(row.names(rna),row.names(auc)),]
  auc <- auc[intersect(row.names(rna),row.names(auc)),]
  print(dim(rna))
  print(dim(auc))
  descriptor <- "rf" #method #"coxrf"
  dump_features(rna, path = paste0("features/",directory,"/",descriptor,"_features.RData"))
  dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response.RData"))
  
  auc_pre <- auc
  for(args in 1:numberofargs){
    if(numberofargs ==1){
      auc <- auc_pre
    } else {
      auc <- cut_df(auc_pre, numberofargs,args)
    }
    dump_features(auc, path = paste0("features/",directory,"/",descriptor,"_response_",as.character(args),".RData"))
  }
}