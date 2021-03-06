setwd("D:/Fangshan/R scripts")
source("Fetch SNPs.R")
library(dplyr)
library(readxl)
library(rio)
library(installr)
Sys.setenv("R_ZIPCMD" = "C:\\Rtools\\bin\\zip.exe")
#test
snps <- c("rs1778155","rs10489177")
Get.SNPinfo(snps,"ALL")
#snp2013
SNPMAF.2013 <- Get.SNPinfo(snp.name.2013,"CHB")

#my snps
my.snps <- read_excel(path="D:\\Fangshan\\SNPS.xlsx", col_names = FALSE)[[1]]
my.snps.info <- Get.SNPinfo(my.snps,"CHB")
my.snps.info$MAF <- as.numeric(pmin(my.snps.info$allele1_freq,my.snps.info$allele2_freq,na.rm=TRUE)) %>% round(digits=3)
#year
my.snps.info$test_2013 <- as.integer(my.snps.info$snpid %in% snp.name.2013)
my.snps.info$test_2016 <- as.integer(my.snps.info$snpid %in% snp.name.2016)
my.snps.info$test_2016pre <- as.integer(my.snps.info$snpid %in% snp.name.2016.pre)
#write to excel table
export(my.snps.info,file = "D:\\Fangshan\\�ҵ�SNP��Ϣ.xlsx", overwrite=TRUE)


