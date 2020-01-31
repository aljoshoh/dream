library(dplyr)

AnvSigGen = function(rnaseq, auc){
  
  #rnaseq: matrix of gene expression data, with samples as rows and genes as columns
  #auc: matrix of drug data, samples as rows and drug as columns
  
  if (any(rownames(rnaseq) != rownames(auc))) {
    auc = auc[rownames(rnaseq), ]
  }
  
  top_expr = colSums(rnaseq) %>% sort(decreasing = T) %>% head(n= 5000)
  rnaseq = rnaseq[, names(top_expr)]
  
  gen_sig = list()
  
  for (d in colnames(auc)) {
    
    auc_d = auc[,d]  %>% na.omit()
    
    if (length(auc_d) < 50){
      next
    }
    auc_d = auc_d[order(auc_d, decreasing = T)]
    top_d = head(auc_d, n=20)  
    bot_d = tail(auc_d, n=20) 
    p_val_vec = vector(mode = "numeric")
    
    for (g in colnames(rnaseq)) {
      rnaseq_g = rnaseq[, g] 
      names(rnaseq_g) = rownames(rnaseq)
      tmp_df = data.frame(expr = c(rnaseq_g[names(top_d)], rnaseq_g[names(bot_d)]),
                          group = rep(c("top", "bot"), each = 20))
      anv = aov(expr ~ group, tmp_df)
      p_val = summary(anv)[[1]][1,5]
      p_val_vec = c(p_val_vec,p_val)

    }
    
    p_val_vec = p.adjust(p_val_vec,"BH")
    index = which(p_val_vec < 0.05)
    gen_sig_d = colnames(rnaseq)[index]
    
    if (length(gen_sig_d) > 0) {
      gen_sig[[d]] = gen_sig_d
    }
  }
  
  sig = selectGen(gen_sig)
  rnaseq_sub = rnaseq[, sig] 
  
  return(rnaseq_sub)
  
}

selectGen = function(gen_sig_list) {
  
  #gen_sig_list: a list of significant genes that explain AUC variation for each drug
  
  gen_selected = vector(mode = "character")
  gen_sig_vec = unlist(gen_sig_list)
  for (i in 1:length(gen_sig_vec)){
    gen = gen_sig_vec[i]
    if (length(gen_sig_vec[gen_sig_vec == gen]) >= 10) {
      gen_selected = c(gen_selected, gen)
    }
  }
  gen_selected = unique(gen_selected)
  return(gen_selected)
}



