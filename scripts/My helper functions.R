# Name: helper functions
# Author: Xiao Xiang
# E-mail:doublex1990@live.cn
# Created Date: 2017-05-09
# Modified Data: 2017-05-09

require(stringr)
require(dplyr)
require(purrr)
library(labelled)
library(rlang)

ColumnToVector <- function(pasted, simplify = TRUE){
  #Convert pasted table into R vectors
  #pasted: data copied from tables
  #simplify: if TRUE, returns a character matrix, otherwise returns a list of vectors
  #value: A list of character Vectors
  pasted.vec <- pasted %>% str_trim() %>% str_split('\n', simplify = TRUE) %>% 
    str_split('\t', simplify = simplify)
  return(pasted.vec)
}

NumToVec <- function(num,len){
  # convert number to index vector
  # num: number to convert
  # len: length of vector
  # value: a vector
  vec <- vector("numeric",len)
  vec[num] <- 1
  return(vec)
}

IdToMat <- function(id){
  #convert ordered id into matrix M
  # id: id number to convert
  # value: a matrix
  ordered.id <- sort(id)
  new.id <- match(ordered.id,id)
  #empty matrix
  mat <- matrix(sapply(new.id,NumToVec, len = length(id)),nrow = length(id))
  return(mat)
}

CharToSciNum <- function(char){
  # convert scientific notation character into numerics
  # char: character
  # value: a numerical vector
  SciNum <- function(x) {
    x <- EmptyRm(x)
    if (length(x) == 2){
      a = x[1] %>% as.numeric()
      base = str_split(x[2],'-',simplify = TRUE)
      b = base[1] %>% as.numeric()
      c = base[2] %>% as.numeric()
      return(a*b^(-c))
    }
    else{
      x <- as.numeric(x)
      return(x[1]*x[2]^x[3])
    }
  }
  if (is.numeric(char)) {
    return(char)
  }
  char.trans <- str_split(char,'[^0-9-.]') %>% sapply(SciNum)
  return(char.trans)
}

EmptyRm <- function(vec){
  # remove empty elements from charactor vector
  # vec: input vector
  # value: a character vector
  return(vec[vec != ""])
} 

## a function to peform two dataframe join
TwoDFMerge <- function(x,y, by = key){
  ## x and y are data.frame objects to merge
  ## by is the key variable
  ## return a data.frame object
  upDate <- function(x.var,y.var){
    df <- data.frame(x.var=x.var,y.var=y.var, stringsAsFactors = FALSE) 
    df <- mutate(df, update=if_else(is.na(y.var), x.var, y.var))
    return(df$update)
  }
  # common and distinct variables
  commonVariable = intersect(colnames(x),colnames(y)) %>% 
    setdiff(by)
  distinct.x = colnames(x) %>% setdiff(commonVariable)
  distinct.y = colnames(y) %>% setdiff(commonVariable)
  # update common variables
  updated <- map2_df(x[commonVariable],y[commonVariable],upDate)
  data <- full_join(x[distinct.x],y[distinct.y], by = by)
  data <- data %>% bind_cols(updated)
  return(data)
}

PriorityMerge <- function(..., key, priority = "default"){
  ## merge dataset with priority, higher priority replace lower priorty
  ## ...: a list of dataframes
  ## priority: the priority of the datasets, lowest to highest by default
  ## value: a dataframe
  require(dplyr)
  args <- list(...)
  n <- length(args) #number of the dataframes
  if (n == 1){
    message("Only one dataframe is specified, return the dataframe without merging")
    return(args[1])
  }
  if (priority == "default"){
    priority <- 1:n 
  }
  args <- args[priority]
  merged.data <- reduce(args, TwoDFMerge, by = key)
  return(merged.data)
}

## generate data dictionary for all categorical variables stored in characters
dicGenerate <- function(data, tolerance=0.1, exclude=NULL){
  ## function that return a tibble containing all levels
  labeller <- function(column, data){
    label <- data[[column]] %>% 
      str_trim() %>% 
      table() %>% 
      as_tibble() %>% 
      select(1) 
    names(label) <- "labels"
    label <- arrange(label, labels)
    if (!is.null(exclude)){
      exclude <- paste0("date|day|time|address",exclude,sep = "|")
    } else {
      exclude <- "date|day|time|address|sign|rs"
    }
    if (nrow(label) >= tolerance*length(na.omit(data[[column]])) |
        str_detect(column,exclude)){
      return(NULL)
    } else {
      label$var_name <- column
      return(label[,c(2,1)])
    }
  }
  character.index <- map_lgl(data, is.character)
  character.variables <- colnames(data)[character.index]
  dict <- map_dfr(character.variables, labeller, data=data)
  return(dict)
}

## change character value to value labels in a dataframe
Value2Label <- function(data,dict){
  vars <- unique(dict$var_name)
  transform <- function(var, data){
    from <- getVariable(dict, labels, var_name == !!var)
    to <- getVariable(dict, value, var_name == !!var)
    named_vec <- setNames(to, from)
    trans <- data[[var]] %>% plyr::mapvalues(from, to) %>% as.numeric()
    trans <- labelled(trans, named_vec)
    return(trans)
  }
  data[vars] <- map(vars, transform, data = data)
  return(data)
}

## get variable from dataframe based on condition
getVariable <- function(data, colname, expression=TRUE){
  colname <- enquo(colname)
  expr <- enquo(expression)
  # If argument is not string, capture and transform it.
  if (!is.character(colname)){
    colname <- quo_name(colname)
  }
  # return selected column in vector
  vector <- dplyr::filter(data, UQE(expr)) %>% `[[`(colname)
  return(vector)
}

## load a object only if not exist
loadObj <- function(obj_name, path,...){
  require(rio, quietly = TRUE)
  obj_name <- enquo(obj_name)
  if (!exists(deparse(UQE(obj_name)))){
    eval(expr(UQE(obj_name) <<- import(path, ...)))
  }
}



