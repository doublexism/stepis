## sourcing external scripts
source("scripts\\external scripts\\HWEtesterror.R")
## functions
library(genetics)
## get AlleleFrequency from a column of snps
getAlleleFrequency <- function(SNP,data,drop.freq=TRUE){
  Alleles <- data[[SNP]]
  Alleles <- str_split(Alleles,':') %>% unlist()
  frequency <- table(Alleles) %>% as.tibble()
  if (nrow(frequency) == 0) {
    return(tibble(rsid = SNP))
  } 
  if (nrow(frequency) == 1){
    a <- tibble(Alleles = 'na', n = 0)
    frequency <- bind_rows(frequency,a)
  }
  frequency <- frequency %>% mutate(AlleleFrequency = round(n/sum(n),3)) 
  if (drop.freq){
    frequency <- frequency %>% select(-n) %>% spread(Alleles, AlleleFrequency)
  } else {
    freq_table <- frequency[c("Alleles", "AlleleFrequency")] %>% spread(Alleles, AlleleFrequency)
    n_table <- frequency[c("Alleles", "n")] %>% mutate(Alleles = paste0(Alleles,"_N")) %>% 
      spread(Alleles, n)
    frequency <- bind_cols(freq_table,n_table)
  }
  if (str_detect(SNP, "_")){
    frequency$rsid <- str_split(SNP, '_',simplify = TRUE)[1]
  } else {
    frequency$rsid <- SNP
  }
  N <- ncol(frequency)
  return(frequency[c(N,seq(1,N-1))])
}

## get genotype frequency
genoFreq <- function(SNP, data){
  genotypes <- unique(data[[SNP]])
  geno_freq <- table(data[[SNP]]) %>% as_tibble()
  N <- nrow(data)
  if (nrow(geno_freq) != 0){
  colnames(geno_freq) <- c("genotype","Count")
  } else {
    geno_freq <- tibble(genotype=NA,Count = 0)
  }
  geno_freq <- geno_freq %>% mutate(rsid = SNP, miss_freq = round((N - sum(Count))/N,3))
  geno_freq <- geno_freq %>% select(rsid,genotype,miss_freq)
  return(geno_freq)
}

## unify the format of the alleles
unifyAlleleFormmat <- function(x){
  unify <- function(element){
    if (is.na(element)){
      return(NA)
    } else if (str_length(element)==1){
      return(paste0(element,':',element))
    } else if (str_length(element==2)){
      return(paste0(str_sub(element,1,1),":",str_sub(element,2,2)))
    }
  }
  if (sum(str_detect(x,':'), na.rm = TRUE) == 0){
    x <- map_chr(x, unify)
  }
  return(x)
}

## unify zygotes
unifyZygotes <- function(x){
  sortPaste <- function(vector){
    vector <- sort(vector)
    vector <- paste0(vector, collapse = ':')
      if (vector == ""){
        return(NA)
      }
    return(vector)
  }
  Alleles <- str_split(x,':')
  x <- map_chr(Alleles, sortPaste)
  return(x)
}
  
## get SNP variables in a dataframe
getSNPName <- function(data){
  var_name <- colnames(data)
  SNP_name <- str_subset(var_name, 'rs[0-9]+')
  return(SNP_name)
}
## Alleles
getAllele <- function(SNP){
  allele <- str_split(SNP,':') %>% unlist() %>% na.omit() %>% unique()
  return(allele)
} 

## get genotype matrix
geno_mat <- function(SNP){
  Allele <- getAllele(SNP)
  return(str_count(SNP,Allele[1]))
}

## add parents id
addPid <- function(data, pid, famid, memid){
  famid <- enquo(famid)
  pid <- enquo(pid)
  memid <- enquo(memid)
  # define a summrize function
  getId <- function(pid, memid, number){
    set.seed(20171114)
    id <- pid[memid  == number]
    proband <- pid[memid == 1]
    if (length(id) > 0){
      return(id)
    } else {
      id_gen <- paste0(proband,"_",number)
      return(id_gen)
    }
  }
  data <- data %>% group_by(!!famid) %>% 
    mutate(motherid = ifelse(memid %in% c(2,3), NA, getId(!!pid, !!memid, 3)),
              fatherid = ifelse(memid %in% c(2,3), NA, getId(!!pid, !!memid, 2)))
  return(data)
}

## summarize genotypes
matchTab <- function(data, famid, SNP){
  famid <- sym(famid)
  SNP <- sym(SNP)
  dat <- data %>% 
    group_by(!!famid) %>% 
    mutate(nmiss = sum(is.na(!!SNP))) %>%
    filter(nmiss == 0) %>%
    summarise(pattern = paste0(UQ(SNP), collapse = ","))
  # explicitly define rownames
  names <- outer(0:2, 0:2, paste, sep = ",") %>% c()
  sum <- table(dat$pattern) %>%
    `[`(names) %>%
    c() %>% 
    matrix(3,3)
  sum[is.na(sum)] <- 0
  return(sum)
}

## hwe test
hwe_test <- function(SNP, data, famid, y = 0.01, u = 0.01){
  num <- matchTab(data, famid, SNP)
  test <- HWEtest.error(num, y, u)
  hwe <- tibble(rsid = SNP, chi = test[[1]], p = test[[2]])
  return(hwe)
}
## hwe test for the unrelated
hwe_unrelate <- function(SNP, data){
  df <- data.frame(chi=NA,pvalue=NA)
  try({g <- genotype(data[[SNP]],sep = ":")
       test <- HWE.test(g, exact =FALSE)
       p <- test$test$p.value
       chisq <- test$test$statistic
       df <- data.frame(chi=chisq,pvalue=p)}, silent = TRUE)
  return(df)
}


## translate genotype to sage format
sageSNP <- function(SNP){
  SNP <- str_replace_all(SNP, ':', "/")
  SNP[is.na(SNP)] <- "./."
  return(SNP)
}
## translate other variables to sage format
sageVar <- function(var){
  var[is.na(var)] <- '.'
  return(var)
}
