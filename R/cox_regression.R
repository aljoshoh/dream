############################################################################################################
# FUNCTIONS
############################################################################################################


library(dplyr)
library(survival)
library(cometExactTest)


plot_cox <- function(
  ### Plots the curves and gives resume statistics
  ### To plot: autoplot(survfit(cox$1,2,3,4))
  ### TODO 
  ### 
  ######################################################
  features = NULL, # feature vectors, $bem_cox
  response = NULL, # response dataframe with survival and covariates, $sav_cox
  feature_name = NULL, # like KRAS, a molecular marker
  strat = NULL, # what feature to filter for in response
  survival = "OS", # also PFS
  covariates = c("1"), # array of covariates
  N = N # minimal number of cell lines per condition
){
  if (all(row.names(features) == row.names(response)) == F) {
    stop("Rownames do not match !")
  }
  
  # Filter for strat variable
  strats <- response[,strat]
  strats_levels <- levels(factor(strats))
  howmany_strats <- length(strats_levels)
  if (nrow(response)==0) {
    stop("Stratification for stat variable failed, not present in dataset ?")
  }
  
  List <- list()
  covariates <- paste("+", paste(covariates, collapse = "+"),sep="")
  for(level in strats_levels){
    logical <- strats == level
    logical[is.na(logical)] <- FALSE # remove NA observations from strat observations
    temp_response <- response[logical,]
    temp_features <- features[logical,]
    if (all(row.names(temp_features) == row.names(temp_response)) == F) {
      stop("Rownames do not match in between strat levels!")
    }
    temp_features <- temp_features[,feature_name]
    temp_response$biomarker <- temp_features
    n1 <- 2#n1 <- length(which(temp_features == "1"))
    n2 <- 2#n2 <- length(which(temp_features == "0"))
    temp_response$covariates <- rep(1, length(temp_response$biomarker))
    if(!(n1 < N | n2 < N)){ # if not minimal amount of mutants are present in one of the investigated populations
      if(survival == "overallSurvival"){
        cox <- coxph(as.formula(paste("Surv(overallSurvival, vitalStatus) ~ biomarker",covariates, sep="")), data = temp_response)
        km_fit <- survfit(as.formula(paste("Surv(overallSurvival, vitalStatus) ~ biomarker","", sep="")), data = temp_response) # add covariates here in the ""
      }
      if(survival == "PFS"){
        cox <- coxph(as.formula(paste("Surv(PFS.months, PFS.event) ~ biomarker",covariates,sep="")), data = temp_response)
        km_fit <- survfit(as.formula(paste("Surv(PFS.months, PFS.event) ~ biomarker","",sep="")), data = temp_response)
      }
      if(survival != "PFS" & survival != "overallSurvival"){
        cox <- NULL; km_fit <- NULL}
    }else{
      cox <- NULL; km_fit <- NULL
    }
    
    List[[level]] <- list(cox=cox, fit=km_fit, strats_levels=level, n1 = n1, n2 = n2)
  }
  
  return(List)
}######################################################




calc_asso_cox <- function(
  ### makes association tests
  ### plots volcano plots from whole dataframe
  ### 
  ### TODO 
  ### 
  ######################################################
  plot_volcano = F,
  features = NULL, # feature vectors, $bem_cox
  response = NULL, # response dataframe with survival and covariates, $sav_cox
  strat = NULL, # what feature to filter for in response
  treatment = "FOLFIRI+Bevacizumab", # also "FOLFIRI+Cetuximab"
  survival = "OS", # also PFS
  covariates = c("1"), # array of covariates
  save = NULL,
  N = 4 # minimal number of tumors per condition
){
  isFunctionCalled <- F
  for(i in 1:length(colnames(features))){
    print(paste("Start run #",as.character(i),sep=""))
    NAME <- colnames(features)[i]
    
    ### MAKE MODEL
    model <<- plot_cox(features = features, response=response, feature_name = NAME, strat = strat, survival = survival, covariates = covariates, N=N)
    if(isFunctionCalled == F){
      asso_cox <- lapply(1:length(model), function(x) data.frame(matrix(nrow=0,ncol=6)))
      lapply(1:length(model), function(x) colnames(asso_cox[[x]]) <<- c("effectsize","P.Value","cfe","strat","npos","nneg") )
      isFunctionCalled <- T
    }
    print(paste("Fill values for each feature: ",NAME,sep=""))
    for(x in 1:length(model)){
      if(!(model[[x]]$n1 < N | model[[x]]$n2 < N)){
        asso_cox[[x]][i,] <- c(model[[x]]$cox$coefficients["biomarker"],
                               (summary(model[[x]]$cox))$coefficients["biomarker","Pr(>|z|)"],
                               NAME,
                               model[[x]]$strats_levels,
                               model[[x]]$n1,
                               model[[x]]$n2
        )
      }else{
        asso_cox[[x]][i,] <- c(NA,NA,NAME, model[[x]]$strats_levels, model[[x]]$n1, model[[x]]$n2)
      }
    }
    print("-----------------------------------------------")
  }
  lapply(1:length(model), function(x) asso_cox[[x]]$fdr <<- p.adjust(asso_cox[[x]]$P.Value, method = "BH"))
  
  for(x in 1:length(asso_cox)){
    asso_cox[[x]]$gene <- colnames(features)
    asso_cox[[x]]$effectsize <- as.numeric(asso_cox[[x]]$effectsize)
    asso_cox[[x]]$P.Value <- as.numeric(asso_cox[[x]]$P.Value)
    asso_cox[[x]]$fdr <- as.numeric(asso_cox[[x]]$fdr)
    asso_cox[[x]]$npos <- as.numeric(asso_cox[[x]]$npos)
    asso_cox[[x]]$nneg <- as.numeric(asso_cox[[x]]$nneg)
    asso_cox[[x]]$confidence <- unlist(lapply(1:nrow(asso_cox[[x]]), function(y) min(asso_cox[[x]]$npos[y],asso_cox[[x]]$nneg[y]))) # confidence score is the minimum of mutants or wild types
  }
  
  save(asso_cox, file=save)
  
  return(print("Done"))
}######################################################



plot_volcano <- function(
  ### plots volcano plot
  ### TODO
  ### Make it more general and interactive
  ######################################################
  df, 
  Cond, 
  sign, 
  annotations='', 
  pval = 0.05, 
  eff = 5,
  REFERENCE = asso_cox_unstrat[[1]],
  PQ = "P.Value"
){
  
  
  df <- df[Cond,]
  significant <- sign[Cond]
  df_arrow <<- df
  
  ## for each column in asso_cox
  df_arrow$old_effectsize <- NA
  df_arrow$old_P.Value <- NA
  
  for(j in 1:nrow(df_arrow)){
    new <<- df_arrow[j,]
    old <<- REFERENCE[REFERENCE$cfe == new$cfe,c("effectsize",PQ)]
    df_arrow[j,"old_effectsize"] <- old$effectsize
    df_arrow[j,"old_P.Value"] <- old[,PQ]
  }
  print(df_arrow)
  df_arrow <- df_arrow
  df_arrow$old_P.Value <- as.numeric(df_arrow$old_P.Value)
  df_arrow$old_effectsize <- as.numeric(df_arrow$old_effectsize)
  
  if(annotations=='')
    annotations = rep(F, length(sign))
  gene_annotation <- (df[,PQ] <= sort(df[,PQ][significant], decreasing = F)[min(50,length(which(significant)))]) #for greatest significant hits
  p1 <- ggplot(df, aes(x=effectsize, y=-log10(df[,PQ]))) + 
    geom_point(aes(size=df$confidence, color = factor(strat)), shape=16, alpha = 0.5) + #size=df$amount; minus confidence measure (is -log(sd))
    labs(title = paste("Significant association for ","COREAD",sep=""), colour="Stratification", size="Confidence score") +
    xlab("Effect Size") +
    #ylab("-log10(p)") + ylim(c(-0.,10)) + xlim(c(-125,125)) + 
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values=if(length(levels(factor(df_arrow$strat)))==1){"brown"}else{c("1"="orange", "2"="blue", "3"="pink", "4"="darkgreen")}, labels= c("1"="CMS1","2"="CMS2","3"="CMS3","4"="CMS4"))+
    geom_segment(data = df_arrow, aes(x = old_effectsize,y = -log10(old_P.Value),xend = effectsize,yend = -log10(df[,PQ]), color=factor(strat)), arrow=arrow(length = unit(0.01, "npc")))+
    if(length(which(gene_annotation & significant)) !=0){geom_text_repel(aes(label= ifelse( gene_annotation & 
                                                                                              significant,df$gene,'')),force = 1, size = 2)}
  
  return(p1)
  
}######################################################


plot_volcano_2 <- function(
  ### plots volcano plot
  ### TODO
  ### Make it more general and interactive
  ######################################################
  df, Cond, annotations='', annotation.object = NULL, pval = 0.05, eff = 5,
  max_p = 8.5, max_eff = 120, name = "", 
  PQ = "P.Value",
  plot_fdr_cutoff = F
){
  
  if(length(which(sort(na.omit(df$fdr))<0.05)) > 0 ){h0 <- sort(na.omit(df[,PQ]))[max(which(sort(na.omit(df$fdr))<0.05))]}else{h0 <- 999999}
  if(length(which(sort(na.omit(df$fdr))<0.1)) > 0 ){h1 <- sort(na.omit(df[,PQ]))[max(which(sort(na.omit(df$fdr))<0.1))]}else{h1 <- 999999}
  if(length(which(sort(na.omit(df$fdr))<0.2)) > 0 ){h2 <- sort(na.omit(df[,PQ]))[max(which(sort(na.omit(df$fdr))<0.2))]}else{h2 <- 999999}
  
  df <- df[Cond,]
  
  p2 <- p1 + geom_point(data=df, aes(x=effectsize, y=-log10(df[,PQ]), size=df$confidence), shape=16, alpha = 0.5) +
    geom_vline(xintercept = eff, linetype="dashed",color = "black", size=0.3)+
    geom_vline(xintercept = -eff, linetype="dashed",color = "black", size=0.3)+
    #geom_hline(yintercept = -log10(pval), linetype="dashed",color = "black", size=0.3)+
    labs(title = paste("FIRE3",": Significant associations for ",name,sep=""), colour="Stratification", size="Confidence score") +
    xlab("Effect Size") +
    ylab(if(PQ=="P.Value"){"-log10(p)"}else{"-log10(q)"}) + ylim(c(-0.,max_p)) + xlim(c((-1.1*max_eff),1.1*max_eff)) + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  if(plot_fdr_cutoff == T){
    p2 <- p2 + geom_hline(yintercept = -log10(h0), linetype="dashed",color = "black", size=0.3)+
      geom_hline(yintercept = -log10(h1), linetype="dashed",color = "black", size=0.3)+
      geom_hline(yintercept = -log10(h2), linetype="dashed",color = "black", size=0.3)+
      annotate("text", x= -max_eff*0.9, y= -log10(h0)-0.035, label = "FDR=5% ", hjust=1)+
      annotate("text", x= -max_eff*0.9, y= -log10(h1)-0.035, label = "FDR=10%", hjust=1)+
      annotate("text", x= -max_eff*0.9, y= -log10(h2)-0.035, label = "FDR=20%", hjust=1)
  }  
  
  return(p2)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

