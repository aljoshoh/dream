library(dplyr)

setwd("/storage/groups/cbm01/workspace/phong.nguyen_new/dream/")

auc = read.csv("aucs.csv", header = T)
rnaseq = read.csv("rnaseq.csv", header = T, check.names = F)

AnvSigGen = function(rnaseq, auc){
  
  #rnaseq: dataframe of gene expression data, with genes as rows and samples as columns
  #auc: dataframe of drug data, with three columns for sample IDs, drug names and AUC values
  
  if (!("lab_id" %in% colnames(auc)) ||
      !("inhibitor" %in% colnames(auc)) ||
      !("auc" %in% colnames(auc))) {
    do.call(return, list("Missing columns or incorrect column's names for AUC dataframe, 
                         Correct colnames: lab_id, inhibitor and auc"))
  }
  
  if (!("Gene" %in% colnames(rnaseq))) {
    do.call(return, list("Missing column 'Gene' or incorrect colname"))
  }
  
  rownames(rnaseq) = rnaseq$Gene
  
  if ("Symbol" %in% colnames(rnaseq)) {
    rnaseq = rnaseq %>% select(-c("Symbol", "Gene"))
  } else {
    rnaseq = rnaseq %>% select(-"Gene")
  }
  
  top_expr = rowSums(rnaseq) %>% sort(decreasing = T) %>% head(n= 5000)
  rnaseq = rnaseq[names(top_expr),]
  
  gen_sig = list()
  
  for (d in levels(auc$inhibitor)) {
    auc_d = filter(auc, inhibitor == d) %>% droplevels()
    if (nrow(auc_d) < 50){
      next
    }
    auc_d = auc_d[order(auc_d$auc, decreasing = T),]
    top_d = head(auc_d$lab_id, n=20) %>% as.character() 
    bot_d = tail(auc_d$lab_id, n=20) %>% as.character()
    p_val_vec = vector(mode = "numeric")
    
    for (g in rownames(rnaseq)) {
      rnaseq_g = rnaseq[g,] %>% as.numeric()
      names(rnaseq_g) = colnames(rnaseq)
      tmp_df = data.frame(expr = c(rnaseq_g[top_d], rnaseq_g[bot_d]),
                          group = rep(c("top", "bot"), each = 20))
      anv = aov(expr ~ group, tmp_df)
      p_val = summary(anv)[[1]][1,5]
      p_val_vec = c(p_val_vec,p_val)

    }
    
    p_val_vec = p.adjust(p_val_vec,"BH")
    index = which(p_val_vec < 0.05)
    gen_sig_d = rownames(rnaseq)[index]
    
    if (length(gen_sig_d) > 0) {
      gen_sig[[d]] = gen_sig_d
    }
  }
  
  sig = selectGen(gen_sig)
  rnaseq_sub = rnaseq[sig, ] %>% t() %>% as.data.frame()
  
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

test = AnvSigGen(rnaseq, auc)
