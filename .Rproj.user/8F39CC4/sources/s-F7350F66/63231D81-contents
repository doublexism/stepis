allele.check <- "rs1029153
rs2787417
rs3865186
rs10489177
rs161802
rs182052
rs2634074
rs6017317
rs934075
rs1116357
rs11614913
rs1387153
rs243021
rs505922
rs7193343
rs2237892
rs2910164
rs4607103
rs579459
rs67839313
rs12445022"
allele.check <- ColumnToVector(allele.check)[,1]
SNP.2013 <- getVariable(SNP.info,SNPs,batch==2013)
SNP.2017 <- getVariable(SNP.info,SNPs,batch==2017)
SNP.intersect <- intersect(SNP.2013,SNP.2017)
allele.check.intersect <- intersect(allele.check, SNP.intersect)
fsdata.genetic.check <- full_join(fsdata.SNP.2013[c("fsid",allele.check.intersect)],
                                  fsdata.SNP.2017[c("fsid",allele.check.intersect)],
                                  by="fsid") %>%
                        select(1,6,8,13,15) %>%
                        filter(complete.cases(.)) 
colnames(fsdata.genetic.check) <- str_replace(colnames(fsdata.genetic.check),
                                              '(x)(y)', '2013')  
# write csv file
write_csv(fsdata.genetic.check, "audit\\fsdata_genetic_check.csv")

