### Functions wrapping the algorithms, cross-validations and output formatting

suppressPackageStartupMessages({
  library(h2o)
  library(dplyr)
  library(tidyverse)
  library(parallel)
  library(glmnet)
  library(glmnetUtils)
  library(rjson)
  library(randomForest)
  library(caret)
  library(tidyr)
  library(survival)
  library(Hmisc)
  library(randomForestSRC)
  library(survcomp)
})
# library(caret) not yet in image !!!


get_features <- function(
  ### Loads features
  ######################################################
  path = NULL, # RData file of feature matrix, e.g. "features/features.RData"
  matrixfy = T
){
  path <- as.character(path)
  message(paste0('Loading ', path))
  if(matrixfy){
    features <- as.matrix(loadRData(path))
  } else {
    features <- loadRData(path)
  }
  
  return(features)
}######################################################

dump_features <- function(
  ### Dumps features
  ######################################################
  R.object = NULL, # R obeject of feature matrix
  path = NULL # where to save
){
  message(paste0("Saving ",as.character(deparse(quote(var)))," in ",as.character(path))," ...")
  save(R.object, file = as.character(path))

  return("Done exporting object !")
}######################################################

cv <- function(
  ### Create cross-validation based on row names
  ### Creates a list of length 2 with kfold elements, each 
  ######################################################
  feature_matrix = NULL, # Numerical feature matrix
  phenotype_matrix = NULL, # Numerical pheontype matrix 
  kfold = 5, # How many folds
  seed = 123, # Setting a fixed seed
  stack = F
){
  #if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}   #<<<<--------------------------------- if error, check if this still runs without model stacking
  #if(!is.matrix(phenotype_matrix)){stop("The phenotype matrix is not a numerical matrix !")}
  if(!stack){
    if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}
    if(!is.matrix(phenotype_matrix)){stop("The feature matrix is not a numerical matrix !")}
  } else {
    if(!is.list(feature_matrix)){stop("The feature matrix is not a list for model stacking !")}
  }
  message("Check for matched row names...")
  
  if(!is.null(phenotype_matrix) & !stack){
    feature_matrix <- feature_matrix[intersect(row.names(feature_matrix), row.names(phenotype_matrix)),]
    phenotype_matrix <- phenotype_matrix[intersect(row.names(feature_matrix), row.names(phenotype_matrix)),]
  }
  
  if(stack){ # For each of the feature dfs in the stacking list, align the samples:
    aligned <- c()
    for(i in 1:length(feature_matrix)){
      feature_matrix[[i]] <- (feature_matrix[[i]])[intersect(row.names(feature_matrix[[i]]), row.names(phenotype_matrix)),]
      phenotype_matrix <- phenotype_matrix[intersect(row.names(feature_matrix[[i]]), row.names(phenotype_matrix)),]
    }
    feature_matrix <- feature_matrix[[1]] # HACK, assuming that the list of feature-matrices has the same rows for each df in list feature_matrix
  }
  
  message(paste0("Feature matrix is of dimension ",as.character(nrow(feature_matrix)),"x",as.character(ncol(feature_matrix))," !"))
  message(paste0("Phenotype matrix is of dimension ",as.character(nrow(phenotype_matrix)),"x",as.character(ncol(phenotype_matrix))," !"))
  message("Creating feature set names ...")
  
  fold_ids <- cut(1:nrow(feature_matrix), breaks = kfold, labels = F)
  
  size <- floor(((kfold-1)/kfold)*nrow(feature_matrix))
  message(paste0("Fold size is ",as.character(size)))
  message(paste0("Test size is ", as.character(nrow(feature_matrix)-size)))
  
  set.seed(seed = seed)
  feature_matrix <- feature_matrix[sample(nrow(feature_matrix)),]

  test_sets <- lapply(1:kfold, function(x) row.names(feature_matrix[fold_ids == x,,drop = F]))
  train_sets <- lapply(1:kfold, function(x) row.names(feature_matrix[fold_ids != x,,drop = F]))
  
  return(list(test_sets = test_sets, 
              train_sets = train_sets
              ))
}######################################################


make_fit <- function(
  ### Runs cross-validation based on row names
  ### Creates a list of length 2 with kfold elements, each 
  ### TODO
  ###
  ######################################################
  feature_matrix = NULL, # Numerical feature matrix
  phenotype_matrix = NULL, # Numerical pheontype matrix 
  folds = NULL, # fold object created in last step
  FUN = NULL, # Pass filtering function for cv sets internally, default is identity function
  hyperparam = c("alpha"=0.5), # Hyperparamter for methods
  method = "glm", # Method
  cvglm = F, # If glm should be cross-validated for the lambda parameter, e.g. for benchmarking
  seed = 123, # Seed for RF method
  returnFit = T, # Logical if $param should be returned for all folds at fit, otherwise it returns lambda.min
  stack = NULL, # for generally having a unique dataframe for each drug
  args = args
){
  if(!stack){
    if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}
    if(!is.matrix(phenotype_matrix)){stop("The feature matrix is not a numerical matrix !")}
    if(length(folds$test_set) != length(folds$train_set)){stop("Folds of testing and training do not match !")}
  } else {
    if(!is.list(feature_matrix)){stop("The feature matrix is not a list for model stacking !")}
  }
  
  message("Initialize Method ...")
  if(method %in% c("dnn")){
    #port <- 8500 #+as.numeric(args), cannot start at multiple ports
    #h2o.init(port=port)
  }
  
  message("Starting to fit ",method," ...")
  score <- matrix(NA,nrow = length(folds$train_set), ncol = ncol(phenotype_matrix))
  gene_names_filtered <- as.data.frame(matrix(NA,nrow = length(folds$train_set), ncol = ncol(phenotype_matrix)))
  param <- as.data.frame(matrix(NA,nrow = length(folds$train_set), ncol = ncol(phenotype_matrix)))
  
  if(method %in% c("rfsurv","cox")){ # Cox needs special treatment, this puts the data-matrix back to a single column
    survival <- as.data.frame(Surv(event = phenotype_matrix[,1], time = phenotype_matrix[,2]))
    row.names(survival) <- row.names(phenotype_matrix)
    phenotype_matrix <- survival
  }
  
  feature_matrix_prestack <- feature_matrix
  for(i in 1:length(folds$train_set)){ # for all folds
    
    ########PHONG METHOD
    if(!stack){
      message(paste0("Executing feature selection for fold ",as.character(i)))
      test <<- feature_matrix[folds$train_sets[[i]],]
      FILTER_FEATURE_NAMES <- FUN(feature = feature_matrix[feature = folds$train_sets[[i]],], auc = phenotype_matrix[folds$train_sets[[i]],])
      ####################
      message(paste0("-> Length of partial response: ",as.character(length(FILTER_FEATURE_NAMES))," !"))
      if(length(FILTER_FEATURE_NAMES) == 0){stop("NO CORRELATED FEATURES FOUND IN FEATURE SELECTION")}
    }
    
    #########BIG FOR LOOP
    for(j in 1:ncol(phenotype_matrix)){
      
      if(stack){ #if model stacking, get the feature matrix for each drug seperately
        feature_matrix <- feature_matrix_prestack[[j]]
        FILTER_FEATURE_NAMES <- lapply(1:ncol(phenotype_matrix), function(x) colnames(feature_matrix))
      }
      message(paste0("...for fold ",as.character(i)))
      message(paste0("...for drug ", as.character(colnames(phenotype_matrix)[j])))
      message(paste0("...selected ",as.character(length(FILTER_FEATURE_NAMES[[j]]))," features..."))
      
      
      y_name <- as.character(colnames(phenotype_matrix)[j])

      y_train <- (phenotype_matrix[folds$train_sets[[i]],j, drop=F])[!is.na(phenotype_matrix[folds$train_sets[[i]],j]),,drop = F]
      
      if(method %in% c("rfsurv","cox")){ ### HACK-BUGFIX, cause of removing the NA from the phenotype data in the line above
        x_train <- feature_matrix[folds$train_sets[[i]],FILTER_FEATURE_NAMES[[j]]] # add different for each drug
      } else {
        x_train <- feature_matrix[names(y_train[,y_name]),FILTER_FEATURE_NAMES[[j]]] # add different for each drug
      }
      
      x_test <- feature_matrix[folds$test_sets[[i]],FILTER_FEATURE_NAMES[[j]], drop = F]
      y_test <- phenotype_matrix[ folds$test_sets[[i]] ,j]
      
      set.seed(seed)
      #######################METHOD
      if(method == "glm"){
        model <- use_glm(x_train, y_train, x_test, y_test,
                         hyperparam = hyperparam,
                         y_name = y_name,
                         cvglm = cvglm,
                         kfold = length(folds$train_set) # how many internal folds
                         )
      }
      
      if(method == "cox"){
        model <- use_cox(x_train, y_train, x_test, y_test,
                         hyperparam = hyperparam,
                         y_name = y_name,
                         cvglm = cvglm,
                         kfold = length(folds$train_set) # how many internal folds
        )
      }
      
      if(method == "rf"){
        model <- use_rf(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = y_name,
                        seed = seed,
                        kfold = length(folds$train_set) # how many internal folds
                        )
      }
      
      if(method == "dnn"){
        model <- use_dnn(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = as.character(colnames(phenotype_matrix)[j]),
                        seed = seed
        )
      }
      
      if(method == "rfsurv"){
        model <- use_rfSurvival(x_train, y_train, x_test, y_test,
                                hyperparam = hyperparam,
                                y_name = y_name,
                                seed = seed
        )
      }
      
      if(method == "glm_bin"){
        model <- use_glm_binary(x_train, y_train, x_test, y_test,
                         hyperparam = hyperparam,
                         y_name = y_name,
                         cvglm = cvglm,
                         kfold = length(folds$train_set) # how many internal folds
        )
      }
      
      if(method == "rf_bin"){
        model <- use_rf_bin(x_train, y_train, x_test, y_test,
                                hyperparam = hyperparam,
                                y_name = y_name,
                                seed = seed,
                                kfold = length(folds$train_set) # how many internal folds
                                
        )
      }
      ############################
      mse <- mean(model$diff*model$diff)
      
      if(!(method %in% c("rfsurv","cox"))){
        cor <- tryCatch(cor(model$pred, y_test, use = "complete.obs", method = "spearman"), error = function(e){message("no correlation calculated");return(NA)})
        if(is.null(cor)){cor <- NA}
        score[i,j] <- cor
      } 
      if(method %in% c("rfsurv","cox")){
        ground <- as.data.frame(as.matrix(y_test))
        s <- cbind(model$pred, ground)
        cor <- concordance.index(x=s[,1], surv.time=s[,2], surv.event=s[,3])$c.index # concordance index
        print(cor)
        score[[i,j]] <- cor
      }
      if(method %in% c("glm_bin","rf_bin")){
        score[i,j] <- NA
      }
      
      if(returnFit == T){
        param[[i,j]] <- list(model$fit)
      }else{
        param[i,j] <- model$fit$lambda.min
      }
      gene_names_filtered[[i,j]] <- list(FILTER_FEATURE_NAMES[[j]])
    }
    ##################
  }
  
  message("End method ...")
  # gives error since is closes connection for all the workers i guess
  if(method %in% c("dnn")){
    #h2o.shutdown(prompt = F)
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
  }

  return(list(score=score,
              param = param,
              gene_names_filtered = gene_names_filtered,
              cv = folds,
              returnFit = returnFit
              ))
}######################################################



make_fit_whole <- function(
  ### Makes the model fit on the whole data
  ### Creates a list of length 2 with kfold elements, each 
  ### TODO
  ###
  ######################################################
  feature_matrix = NULL, # Numerical feature matrix
  phenotype_matrix = NULL, # Numerical pheontype matrix 
  FUN = function(x,y){res<-lapply(1:ncol(y),function(z) colnames(x));return(res)}, # Pass filtering function for cv sets internally, default is identity function
  hyperparam = c("alpha"=0.5), # Hyperparamter for methods
  method = "glm", # Method
  seed = 123, # Seed for RF method
  stack = F
){
  if(!is.matrix(feature_matrix)){stop("The feature matrix is not a numerical matrix !")}
  if(!is.matrix(phenotype_matrix)){stop("The feature matrix is not a numerical matrix !")}
  if(!stack){
    intersect <- intersect(row.names(feature_matrix), row.names(phenotype_matrix))
  }else{
    intersect <- intersect(row.names(feature_matrix[[1]]), row.names(phenotype_matrix)) #[[1]] just as a test
  }
  
  message("Initialize Method ...")
  if(method %in% c("dnn")){
    h2o.init(port=8502, min_mem_size = "20G")
  }
  
  message("Starting to fit ",method," ...")
  score <- matrix(NA,nrow = 1, ncol = ncol(phenotype_matrix))
  param <- as.data.frame(matrix(NA,nrow = 1, ncol = ncol(phenotype_matrix)))
  gene_names_filtered <- as.data.frame(matrix(NA,nrow = 1, ncol = ncol(phenotype_matrix)))
  
  if(method %in% c("cox","rfsurv")){ # Cox needs special treatment, this puts the data-matrix back to a single column
    survival <- as.data.frame(Surv(event = phenotype_matrix[,1], time = phenotype_matrix[,2]))
    row.names(survival) <- row.names(phenotype_matrix)
    phenotype_matrix <- survival
  }
  
  feature_matrix_prestack <- feature_matrix
  
  ########PHONG METHOD
  #message(paste0("Executing feature selection for overall ",""))
  #test <- feature_matrix[intersect,]
  #FILTER_FEATURE_NAMES <- colnames(FUN(feature_matrix[intersect,], phenotype_matrix[intersect,]))
  ####################
  #message(paste0("-> Found ",as.character(length(FILTER_FEATURE_NAMES))," genes !"))
  
  if(!stack){
    message(paste0("Executing feature selection for overall ",""))
    test <<- feature_matrix[intersect,]
    FILTER_FEATURE_NAMES <- FUN(feature_matrix[intersect,], phenotype_matrix[intersect,])
    ####################
    message(paste0("-> Length of partial response: ",as.character(length(FILTER_FEATURE_NAMES))," !"))
  }
  
  for(j in 1:ncol(phenotype_matrix)){
    message(paste0("...for drug ", as.character(colnames(phenotype_matrix)[j])))
    for(i in 1:1){ # there is only one fold....
      
      if(stack){ #if model stacking, get the feature matrix for each drug seperately
        feature_matrix <- feature_matrix_prestack[[j]]
        FILTER_FEATURE_NAMES <- lapply(1:ncol(phenotype_matrix), function(x) colnames(feature_matrix))
      }
      message(paste0("...selected ",as.character(length(FILTER_FEATURE_NAMES[[j]]))," features..."))
      
      #message(paste0("...for fold ",as.character(i)))
      y_name <- as.character(colnames(phenotype_matrix)[j])
      y_train <- (phenotype_matrix[intersect,j, drop=F])[!is.na(phenotype_matrix[intersect,j]),,drop = F]
      
      if(method %in% c("cox","rfsurv")){ ### HACK-BUGFIX, cause of removing the NA from the phenotype data in the line above
        x_train <- feature_matrix[intersect,FILTER_FEATURE_NAMES[[j]]]
      } else {
        x_train <- feature_matrix[names(y_train[,y_name]),FILTER_FEATURE_NAMES[[j]]]
      }
      
      #x_train <- feature_matrix[names(y_train[,y_name]),]
      
      message(paste0("        Progress ",as.character(round(j/ncol(phenotype_matrix)*100)),"% !"))
      x_test <- NULL #feature_matrix[intersect,]
      y_test <- NULL #phenotype_matrix[ intersect ,j]

      set.seed(seed)
      #######################METHOD
      if(method == "glm"){
        model <- use_glm(x_train, y_train, x_test, y_test,
                            hyperparam = hyperparam,
                            y_name = y_name,
                            cvglm = F
        )
      }
      
      if(method == "cox"){
        model <- use_cox(x_train, y_train, x_test, y_test,
                         hyperparam = hyperparam,
                         y_name = y_name,
                         cvglm = F
        )
      }
      
      if(method == "rf"){
        model <- use_rf(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam[[j]], # !!!! optimal hyperparameters
                        y_name = y_name,
                        seed = seed
        )
      }
      
      if(method == "dnn"){
        model <- use_dnn(x_train, y_train, x_test, y_test,
                        hyperparam = hyperparam,
                        y_name = as.character(colnames(phenotype_matrix)[j]),
                        seed = seed
        )
      }
      
      if(method == "rfsurv"){
        model <- use_rfSurvival(x_train, y_train, x_test, y_test,
                         hyperparam = hyperparam,
                         y_name = y_name,
                         seed = seed
        )
      }
      
      if(method == "glm_bin"){
        model <- use_glm_binary(x_train, y_train, x_test, y_test,
                                hyperparam = hyperparam,
                                y_name = y_name,
                                cvglm = F
        )
      }
      
      if(method == "rf_bin"){
        model <- use_rf_bin(x_train, y_train, x_test, y_test,
                               hyperparam = hyperparam,
                               y_name = y_name,
                               seed = seed
                              
                               
        )
      }
      #######################METHOD
      mse <- NA
      cor <- NA
      
      score[i,j] <- cor
      param[[i,j]] <- list(model$fit) #<- here we can either asses feature weights, or store other things from model$
      gene_names_filtered[[i,j]] <- list(FILTER_FEATURE_NAMES[[j]])
    }
  }
  
  message("End method ...")
  # gives error since is closes connection for all the workers i guess
  if(method %in% c("dnn")){
    #h2o.shutdown(prompt = F)
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
  }
  
  return(list(score=score,
              param = param, # returns the whole fit object
              gene_names_filtered = gene_names_filtered
  ))
}######################################################
