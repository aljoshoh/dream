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



import_cv_results <- function(
  ###
  #####################################################
  partial_path = NULL, # pattern to match the path to the different instances 
  directory = NULL
){
  files <- mixedsort(list.files(pattern = partial_path, path = directory))
  cum_object <- list()
  count <- 1
  for( file in files ){
    object <- loadRData(paste0(directory,"/",file))
    print(count)
    if(count >1){
      cum_object$score <- cbind(cum_object$score, object$score)
      cum_object$param <- cbind(cum_object$param, object$param)
      cum_object$gene_names_filtered <- cbind(cum_object$gene_names_filtered, object$gene_names_filtered)
    } else {
      cum_object <- object
    }
    count <- count +1
  }
  return(cum_object)
}######################################################

























### gtools package function
mixedsort <- function(x,
                      decreasing=FALSE,
                      na.last=TRUE,
                      blank.last=FALSE,
                      numeric.type=c("decimal", "roman"),
                      roman.case=c("upper","lower","both")
)
{
  ord <- mixedorder(x,
                    decreasing=decreasing,
                    na.last=na.last,
                    blank.last=blank.last,
                    numeric.type=numeric.type,
                    roman.case=roman.case
  )
  x[ord]
}
  
mixedorder <- function(x,
                       decreasing=FALSE,
                       na.last=TRUE,
                       blank.last=FALSE,
                       numeric.type=c("decimal", "roman"),
                       roman.case=c("upper","lower","both")
)
{
  # - Split each each character string into an vector of strings and
  #   numbers
  # - Separately rank numbers and strings
  # - Combine orders so that strings follow numbers
  
  numeric.type <- match.arg(numeric.type)
  roman.case   <- match.arg(roman.case)
  
  if(length(x)<1)
    return(NULL)
  else if(length(x)==1)
    return(1)
  
  if( !is.character(x) )
    return( order(x, decreasing=decreasing, na.last=na.last) )
  
  delim="\\$\\@\\$"
  
  if(numeric.type=="decimal")
  {
    regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"  # uses PERL syntax
    numeric <- function(x) as.numeric(x)
  }
  else if (numeric.type=="roman")
  {
    regex <- switch(roman.case,
                    "both"  = "([IVXCLDMivxcldm]+)",
                    "upper" = "([IVXCLDM]+)",
                    "lower" = "([ivxcldm]+)"
    )
    numeric <- function(x) roman2int(x)
  }
  else
    stop("Unknown value for numeric.type: ", numeric.type)
  
  nonnumeric <- function(x)
  {
    ifelse(is.na(numeric(x)), toupper(x), NA)
  }
  
  x <- as.character(x)
  
  which.nas <- which(is.na(x))
  which.blanks <- which(x=="")
  
  ####
  # - Convert each character string into an vector containing single
  #   character and  numeric values.
  ####
  
  # find and mark numbers in the form of +1.23e+45.67
  delimited <- gsub(regex,
                    paste(delim,"\\1",delim,sep=""),
                    x,
                    perl=TRUE)
  
  # separate out numbers
  step1 <- strsplit(delimited, delim)
  
  # remove empty elements
  step1 <- lapply( step1, function(x) x[x>""] )
  
  # create numeric version of data
  suppressWarnings( step1.numeric <-  lapply( step1, numeric ) )
  
  # create non-numeric version of data
  suppressWarnings( step1.character <- lapply( step1, nonnumeric ) )
  
  # now transpose so that 1st vector contains 1st element from each
  # original string
  maxelem <- max(sapply(step1, length))
  
  step1.numeric.t <- lapply(1:maxelem,
                            function(i)
                              sapply(step1.numeric,
                                     function(x)x[i])
  )
  
  step1.character.t <- lapply(1:maxelem,
                              function(i)
                                sapply(step1.character,
                                       function(x)x[i])
  )
  
  # now order them
  rank.numeric   <- sapply(step1.numeric.t, rank)
  rank.character <- sapply(step1.character.t,
                           function(x) as.numeric(factor(x)))
  
  # and merge
  rank.numeric[!is.na(rank.character)] <- 0  # mask off string values
  
  rank.character <- t(
    t(rank.character) +
      apply(matrix(rank.numeric),2,max,na.rm=TRUE)
  )
  
  rank.overall <- ifelse(is.na(rank.character),rank.numeric,rank.character)
  
  order.frame <- as.data.frame(rank.overall)
  if(length(which.nas) > 0)
    if(is.na(na.last))
      order.frame[which.nas,] <- NA
  else if(na.last)
    order.frame[which.nas,] <- Inf
  else
    order.frame[which.nas,] <- -Inf
  
  if(length(which.blanks) > 0)
    if(is.na(blank.last))
      order.frame[which.blanks,] <- NA
  else if(blank.last)
    order.frame[which.blanks,] <- 1e99
  else
    order.frame[which.blanks,] <- -1e99
  
  order.frame <- as.list(order.frame)
  order.frame$decreasing <- decreasing
  order.frame$na.last <- NA
  
  retval <- do.call("order", order.frame)
  
  return(retval)
}
  

  