
AnvSigNumFeature = function(feature, auc){
  
  #feature: matrix of numerical data, with samples as rows and features as columns
  #auc: matrix of drug data, samples as rows and drug as columns
  
  require(dplyr)
  require(tidyr)
  
  if (any(rownames(feature) != rownames(auc))) {
    auc = auc[rownames(feature), ]
  }
  
  feature_sig = list()
  
  for (d in colnames(auc)) {
    
    auc_d = auc[,d]  %>% na.omit()
    
    if (length(auc_d) < 20){
      next
    }
    auc_d = auc_d[order(auc_d, decreasing = T)]
    top_d = head(auc_d, n=10)  
    bot_d = tail(auc_d, n=10) 
    p_val_vec = vector(mode = "numeric", length = ncol(feature))
    names(p_val_vec) = colnames(feature)
    
    for (f in colnames(feature)) {
      feature_f = feature[, f] 
      names(feature_f) = rownames(feature)
      tmp_df = data.frame(value = c(feature_f[names(top_d)], feature_f[names(bot_d)]),
                          group = rep(c("top", "bot"), each = 10))
      anv = aov(value ~ group, tmp_df)
      p_val = summary(anv)[[1]][1,5]
      p_val_vec = c(p_val_vec,p_val)
      
    }
    
    p_val_vec = p.adjust(p_val_vec,"BH")
    p_val_vec = p_val_vec[order(p_val_vec)]
    index = names(head(p_val_vec, n=100))
    feature_sig_d = feature[,index]
    
    feature_sig[[d]] = feature_sig_d
  }
  
  return(feature_sig)
  
}

AnvSigCatFeature = function(feature, auc){
  
  #feature: matrix of categorical feature data, with samples as rows and features as columns
  #auc: matrix of drug data, samples as rows and drug as columns
  
  require(dplyr)
  require(tidyr)
  
  if (any(rownames(feature) != rownames(auc))) {
    auc = auc[rownames(feature), ]
  }
  
  feature_sig = list()
  
  for (d in colnames(auc)) {
    
    auc_d = auc[,d]  %>% na.omit()
    
    if (length(auc_d) < 20){
      next
    }
     
    p_val_vec = vector(mode = "numeric", length = ncol(feature))
    names(p_val_vec) = colnames(feature)
    
    for (f in colnames(feature)) {
      feature_f = feature[, f] 
      names(feature_f) = rownames(feature)
      feature_f = feature_f[names(auc_d)]
      
      if (nlevels(as.factor(feature_f)) < 2) {
        next
      }
      
      tmp_df = data.frame(auc = auc_d,
                          group = feature_f)
      anv = aov(auc ~ group, tmp_df)
      p_val = summary(anv)[[1]][1,5]
      p_val_vec = c(p_val_vec,p_val)
      
    }
    
    p_val_vec = p.adjust(p_val_vec,"BH")
    p_val_vec = p_val_vec[order(p_val_vec)]
    index = names(head(p_val_vec, n=100))
    feature_sig_d = feature[,index]
    
    feature_sig[[d]] = feature_sig_d
    
  }
  
  return(feature_sig)
  
}


gen_select = AnvSigNumFeature(feature = rnaseq, auc = auc)

  