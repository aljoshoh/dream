
AnvSigNumFeature = function(feature, auc){
  
  #feature: matrix of numerical data, with samples as rows and features as columns
  #auc: matrix of drug data, samples as rows and drug as columns
  
  require(dplyr)
  require(tidyr)
  
  if (any(rownames(feature) != rownames(auc))) {
    auc = auc[rownames(feature), ]
  }
  
  top = apply(feature, 2, sd) %>% sort(decreasing = T) %>% head(n=5000)
  feature = feature[, names(top)]
  
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
      feature_f = feature_f[names(auc_d)]
      
      mod = lm(auc_d ~ feature_f)
      
      if (dim(coef(summary(mod))) == c(2,4)){
        p_val = coef(summary(mod))[2,4]
        p_val_vec[f] = p_val
      } else {
        p_val_vec[f] = NA
      }
      
    }
    
    p_val_vec = p.adjust(p_val_vec,"BH")
    p_val_vec = p_val_vec[order(p_val_vec)]
    
    if (length(p_val_vec) > 100){
      feature_sig[[d]] = names(head(p_val_vec, n=100))
    } else {
      feature_sig[[d]] = names(p_val_vec) # fixed a bug here
    }
    
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
      
      if (nlevels(as.factor(feature_f)) < 2 |
          any(table(as.factor(feature_f)) < 1)) {
        next
      }
      
      tmp_df = data.frame(auc = auc_d,
                          group = feature_f)
      anv = aov(auc ~ group, tmp_df)
      p_val = summary(anv)[[1]][1,5]
      p_val_vec[f] = p_val
      
    }
    
    p_val_vec = p.adjust(p_val_vec,"BH")
    p_val_vec = p_val_vec[order(p_val_vec)]
    
    feature_sig[[d]] = names(head(p_val_vec, n=5))
    
  }
  
  return(feature_sig)
  
}

  
