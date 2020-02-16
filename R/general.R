cut_df <- function(
  ### cuts a dataframe in columns and returns the one with a certain index
  ######################################################  
  df = NULL,
  totaln = NULL,
  index = NULL
){
  tmp <- cut(1:ncol(df), breaks = totaln, labels = F)
  df <- df[,tmp == index]
  return(df)
}######################################################


loadRData <- function(
  ### Loads an RData file, and returns it
  ######################################################
  fileName
){
  load(fileName)
  get(ls()[ls() != "fileName"])
}######################################################


get_params <- function(
  ### Loads hyperparameters
  ######################################################
  json.file = NULL # json file of paramters, e.g. "params/params.json"
){
  json.file <- as.character(json.file)
  message(paste0('Loading ', json.file))
  params <- fromJSON(json.file)
  
  return(params)
}######################################################


combine_results <- function(
  ### Loads results and saves the full df
  ######################################################
  name = NULL # path until the instance variable
){
  
  
}
  
  
  
  
  
  
  