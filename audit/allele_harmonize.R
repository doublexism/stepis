### read the alleles to be harmonized
allele_harmonize <- "rs11614913	k
rs161802	k
rs340874	r
rs505922	k
rs6017317	k
rs1029153	kr
rs1116357	kr
rs1387153	k
rs2107595	r
rs2237892	k
rs2240419	r
rs243021	kr
rs2540489	r
rs2634074\tk
rs2787417	k
rs2910164	k
rs3865186	kr
rs4607103	k
rs4812829	d
rs501120	r
rs556621	r
rs67839313	k
rs7193343	k
rs12445022	k
rs182052	k
rs579459	k
"
fsdata.SNP.2017 <- fsdata[c("fsid","pid", SNP.names.2017)]
colnames(fsdata.SNP.2017) <- c("fsid","pid",getVariable(SNP.info, SNPs, variable %in% SNP.names.2017))
allele_harmonize <- ColumnToVector(allele_harmonize, simplify = FALSE)
allele_names <- unlist(allele_harmonize)[seq(1,length(allele_harmonize)*2,2)]
snpHarmonizer <- function(param){
  rsid <- param[1]
  action <- param[2]
  print(rsid)
  print(action)
  vector <- fsdata.SNP.2017[[rsid]]
  SNP_mat <- str_split(vector,':') %>% unlist()%>% na.omit() %>% unique()  
  dict <- tibble(from = c("A","T","G","C"))
  if (str_detect(action,'k')){
    old_SNP <- paste0(SNP_mat, collapse = "")
    new_SNP <- paste0(SNP_mat[2:1], collapse = "")
   
    vector <- chartr(old_SNP, new_SNP, vector)
  }
  if (str_detect(action,'r')){
    vector <- chartr("AGTC","TCAG", vector)
  }
  if (action == 'd'){
    return(NULL)
  }
  return(vector)
}

## harmonize the 2017 SNPs
fsdata.SNP.2017[allele_names] <- map(allele_harmonize, snpHarmonizer)

