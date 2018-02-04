## read in the raw data
fsdata.before.2017 <- read_csv(file = "data\\fsdata_170219.csv",
                               locale = locale(encoding='GB2312'),
                               na = c('NA','?','','.','missing'),
                               guess_max = 10000)
fsdata.genotype.2017 <- read_csv(file = "data\\fsdata_170721.csv",
                                 locale = locale(encoding='GB2312'),
                                 na = c('NA','?','','.','missing'),
                                 guess_max = 10000)
fsdata.genotype.2017 <- fsdata.genotype.2017 %>% filter(!is.na(fsid))

# unify variant call
SNP.names.before.2017 <- str_subset(colnames(fsdata.before.2017),"rs[0-9]+_[0-9](pre)?")
fsdata.before.2017[SNP.names.before.2017] <- map(fsdata.before.2017[SNP.names.before.2017]
                                                 ,unifyAlleleFormmat)

## check genotyping errors
SNPs <- colnames(fsdata.genotype.2017[2:ncol(fsdata.genotype.2017)])
allele.frequency <- map_dfr(SNPs, getAlleleFrequency, data=fsdata.genotype.2017)
# get 1kg genotype info from website
loadObj(allele.1kg, "data\\allele_1kg.rds")

rs.id <- str_split(SNPs,'_',simplify = TRUE)[,1]
allele.1kg <- Get.SNPinfo(rs.id,pop="CHB")
write_rds(allele.1kg, "data\\allele_1kg.rds")
allele.check <- left_join(allele.frequency, allele.1kg,by=c("rsid"="SNP"))
write_csv(allele.check, "audit\\allele_check.csv")

## merge datasets
fsdata <- left_join(fsdata.before.2017, fsdata.genotype.2017,by = "fsid")

## variables
varname <- colnames(fsdata)
SNP.names <- str_subset(varname,"rs[0-9]+_[0-9](pre)?")
SNP.info <- data.frame(variable = SNP.names, stringsAsFactors = FALSE)
SNP.info$SNPs <- str_split(SNP.names, '_', simplify = TRUE)[,1]
SNP.info$batch <- str_split(SNP.names, '_', simplify = TRUE)[,2]
SNP.names.2013 <- getVariable(SNP.info, variable, batch=="2013")
SNP.names.2016 <- getVariable(SNP.info, variable, batch=="2016")
SNP.names.2016pre <- getVariable(SNP.info, variable, batch=="2016pre")
SNP.names.2017 <- getVariable(SNP.info, variable, batch=="2017")

## get variable name and type
fsdata.dict <- data.frame(names=varname,
                          types=map_chr(fsdata,class),
                          row.names = NULL)

## merge genotyping data
fsdata.SNP.2013 <- fsdata[c("fsid","pid",SNP.names.2013)]
## fix rs1801278
fsdata.SNP.2013$rs1801278 <- str_replace_all(fsdata.SNP.2013$rs1801278, "T","A")
fsdata.SNP.2016 <- fsdata[c("fsid","pid",SNP.names.2016)]
fsdata.SNP.2016pre <- fsdata[c("fsid","pid", SNP.names.2016pre)]
fsdata.SNP.2017 <- fsdata[c("fsid","pid", SNP.names.2017)]
colnames(fsdata.SNP.2013) <- c("fsid","pid",getVariable(SNP.info, SNPs, variable %in% SNP.names.2013))
colnames(fsdata.SNP.2016) <- c("fsid","pid",getVariable(SNP.info, SNPs, variable %in% SNP.names.2016))
colnames(fsdata.SNP.2016pre) <- c("fsid","pid",getVariable(SNP.info, SNPs, variable %in% SNP.names.2016pre))
colnames(fsdata.SNP.2017) <- c("fsid","pid",getVariable(SNP.info, SNPs, variable %in% SNP.names.2017))
source("audit\\allele_harmonize.R")
## prioritized merge
fsdata.update <- PriorityMerge(fsdata.SNP.2013, fsdata.SNP.2016pre, fsdata.SNP.2016, fsdata.SNP.2017,
              key = "pid")
## harmonize zygotes
SNP.names.final <- colnames(fsdata.update) %>% setdiff(c("fsid","pid"))

fsdata.update[SNP.names.final] <- map(fsdata.update[SNP.names.final], unifyZygotes)

## check allele frequency again
allele.frequency <- map_dfr(SNP.names.final, 
                            getAlleleFrequency,
                            data=fsdata.update) %>%
  rowwise() %>%
  mutate(Minor_Allele = c("A","T","C","G")[which.min(c(A,T,C,G))],
         MAF = min(c(A,T,C,G), na.rm = TRUE)) %>%
  arrange(rsid)
#output to csv
write_excel_csv(allele.frequency, "audit\\allele_frequency_check.csv")

## remove SNPs in original data and merge new data
SNP.names.all <- str_subset(varname,"rs[0-9]+")
fsdata[SNP.names.all] <- NULL
fsdata.final <- left_join(fsdata, fsdata.update, by="pid")

## label management
if (!exists("fsdata.final")){
  fsdata.final <- read_rds("data\\fsdata_1106.rds")
}
# read the labels into dataset

## Save results
write_rds(fsdata.update,"data\\fsdata.rds")
write_rds(fsdata.final,"data\\fsdata_1201.rds")
write_excel_csv(fsdata.final,"data\\fsdata_1201_UTF8.csv")

fsdata.patch <- fsdata.final %>% select(pid, rs556621, rs2634074)
write_rds(fsdata.patch, "data\\fsdata.patch.rds")
