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