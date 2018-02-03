SNP.2013 <- getVariable(SNP.info, SNPs, batch == "2013")
SNP.2016pre <- getVariable(SNP.info, SNPs, batch == "2016pre")
SNP.2016 <- getVariable(SNP.info, SNPs, batch == "2016")
SNP.2017 <- getVariable(SNP.info, SNPs, batch == "2017")

overlap13 <- intersect(SNP.2013, SNP.2017)
overlap16 <- intersect(SNP.2016pre,SNP.2017)
check <- function(snp,batch){
  var_1 <- paste0(snp, "_",batch)
  var_2 <- paste0(snp, "_2017")
  data <- fsdata[c(var_1,var_2)] %>%
   na.omit()
  if (nrow(data) != 0){
    return(data)
  }
}

data_13 <- map(overlap13, check, batch="2013")
data_16 <- map(overlap16, check, batch="2016pre") %>% Filter(Negate(is.null),.)

table13 <- map(fsdata.SNP.2013[overlap13],table)
table17 <- map(fsdata.SNP.2017[overlap13],table)
print_allele <- function(i){
  print(c(table13[i], table17[i]))
}
