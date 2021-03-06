
AnvSigNumFeature = function(feature, auc){
  
  #feature: matrix of numerical data, with samples as rows and features as columns
  #auc: matrix of drug data, samples as rows and drug as columns
  
  
  if (any(rownames(feature) != rownames(auc))) {
    auc = auc[rownames(feature), ]
  }
  
  top = apply(feature, 2, function(x) sd(x, na.rm = )) %>% sort(decreasing = T) %>% head(n=10000) # alex edit, select 10000
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
      
      #mod = lm(auc_d ~ feature_f)
      mod = cor(auc_d, feature_f, use = "complete.obs")
      
      #if (dim(coef(summary(mod))) == c(2,4)){
      #  p_val = coef(summary(mod))[2,4]
      #  p_val_vec[f] = p_val
      #} else {
      #  p_val_vec[f] = NA
      #}
      p_val_vec[f] = mod
    }
    
    #p_val_vec = p.adjust(p_val_vec,"BH") unnnecesary if we just want the ordering
    p_val_vec = p_val_vec[order(p_val_vec)]
    
    if (length(p_val_vec) > 1000){
      feature_sig[[d]] = names(head(p_val_vec, n=1000)) # alex edit, select 1000
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

  








#### TODO make object for survival !
AnvSigSurvFeature = function(feature, auc, all = F, selection=T,...){ # do not touch the selection parameter
  
  #feature: matrix of numerical data, with samples as rows and features as columns
  #auc: matrix of drug data, samples as rows and drug as columns
  #names(auc) = row.names(feature)
  
  if (any(row.names(feature) != row.names(auc))) {
    auc = auc[row.names(feature), ]
  }
  names(auc) = row.names(feature)
  gex_names <- loadRData("metadata/gex_names.RData")
  feature_gex <- tryCatch(feature[,gex_names], error = function(e) NULL)
  if(!is.null(feature_gex)){
    top = apply(feature_gex, 2, function(x) sd(x, na.rm = T)) %>% sort(decreasing = T) %>% head(n=10000) # alex edit, select 1000
    feature_gex = feature_gex[, names(top)]
  }
  feature <- feature[,!colnames(feature) %in% gex_names]
  if(all == T){feature_gex <- cbind(feature, feature_gex); feature <- NULL}
  
  feature_sig = list()
  
  for (i in 1:1) {
    auc_d = auc 
    
    if (length(auc_d) < 20){
      next
    }
    if(!is.null(ncol(feature_gex)) & selection){
      message("Feature selection is turned on ...")
      p_val_vec = vector(mode = "numeric", length = ncol(feature_gex))
      names(p_val_vec) = colnames(feature_gex)
    
      for (f in colnames(feature_gex)) {
        feature_f = feature_gex[, f]
        feature_f = feature_f[names(auc_d)]
      
        #mod = lm(auc_d ~ feature_f)
        #mod = cor(auc_d, feature_f, use = "complete.obs")
        mod = summary(coxph(auc_d~ feature_f))$coefficients[5]
      
        #if (dim(coef(summary(mod))) == c(2,4)){
          #  p_val = coef(summary(mod))[2,4]
          #  p_val_vec[f] = p_val
        #} else {
        #  p_val_vec[f] = NA
        #}
        p_val_vec[f] = mod
      }
    }else{
      p_val_vec <-1
    }
      #p_val_vec = p.adjust(p_val_vec,"BH") unnnecesary if we just want the ordering
      p_val_vec = p_val_vec[order(p_val_vec)]
    if (length(p_val_vec) > 1000){
      feature_sig[[1]] = c(names(head(p_val_vec, n=1000)),colnames(feature)) # alex edit, select 1000
    } else {
      feature_sig[[1]] = c(names(p_val_vec),colnames(feature)) # fixed a bug here
    }
    
  }
  
  return(feature_sig)
  
}


genesets_survival = function(feature, auc){
set <- loadRData("features/rna-surv/GeneSets_survival.RData")$Genes
return(NULL)
}

