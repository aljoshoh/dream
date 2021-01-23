### Functions to engineer feature sets from binary event matrices (e.g. genomic variants)
#combinations_bem: make all pairwise combinations
#exclusivity_bem: make feature data frame from logical combinations of pairwise combinations


combinations_bem <- function(
  ### names of filter must mach with the names defined here in the function, the sep="_or_"
  ### TODO
  ### think about not combinations
  ### assign permutative significance
  ###
  ######################################################
  bem = bem, # binary event matrix with characters, names do not matter
  logic = "or", # alternatively, "and"
  min_altered = 4,
  filter = NULL # array of names of columnes of $bem
){
  bem <- data.matrix(bem)
  frequency <- apply(bem, 2, function(x) sum(x))
  frequency <- frequency[frequency >= min_altered]
  print(paste("The number of combinations will be: ",as.character(length(frequency)*(length(frequency)-1)/2),sep=""))
  
  bem_min <<- bem[, apply(bem, 2, function(x) sum(x)) > min_altered]
  
  print("Starting pairwise AND/OR calculation...")
  combinations <- combn(1:ncol(bem_min),2)
  nm1 <- unlist(lapply(1:ncol(combinations), function(x) paste(colnames(bem_min)[combinations[1,x]],colnames(bem_min)[combinations[2,x]],sep="_or_")))
  if(logic == "or"){
    res <- lapply(1:ncol(combinations), function(x) unlist(as.character(as.numeric(
      as.logical(bem_min[,combinations[1,x]]) | as.logical(bem_min[,combinations[2,x]]) 
    )))
    )
  }
  if(logic == "and"){
    res <- lapply(1:ncol(combinations), function(x) unlist(as.character(as.numeric(
      as.logical(bem_min[,combinations[1,x]]) & as.logical(bem_min[,combinations[2,x]])
    )))
    )
  }
  res <- as.data.frame(res)
  colnames(res) <- nm1
  row.names(res) <- row.names(bem)
  return(res[,filter])
}######################################################



exclusivity_bem <- function(
  ###
  ### TODO 
  ### think about not combinations
  ### assign permutative significance
  ###
  ######################################################
  bem = bem, # binary event matrix with characters, names do not matter
  logic = "exclusivity", # alternatively, "and"
  min_altered = 4
){
  bem <- data.matrix(bem)
  frequency <- apply(bem, 2, function(x) sum(na.omit(x)))
  frequency <- frequency[frequency >= min_altered]
  print(paste("The number of comparisions will be: ",as.character(length(frequency)*(length(frequency)-1)/2),sep=""))
  
  bem_min <<- bem[, as.logical(apply(bem, 2, function(x) sum(na.omit(x))) > min_altered)]
  
  print("Starting mutual exclusivity calculation...")
  combinations <- combn(1:ncol(bem_min),2)
  nm1 <- unlist(lapply(1:ncol(combinations), function(x) paste(colnames(bem_min)[combinations[1,x]],colnames(bem_min)[combinations[2,x]],sep="_or_")))
  
  if(logic == "exclusivity"){
    res <- lapply(1:ncol(combinations), function(x){
      comet_exact_test(table(bem_min[,combinations[1,x]], bem_min[,combinations[2,x]]) )
    }
    )
  }
  
  res <- unlist(res)
  names(res) <- nm1
  res_adj <- p.adjust(res, method = "BH")
  return(list(pval=res, fdr=res_adj, names=names(res)))
}######################################################