biocLite("AnnotationHub")
biocLite("ensembldb")
library(biomaRt)
library(AnnotationHub)
mart.snp <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp","http://Aug2017.archive.ensembl.org")
listDatasets(mart.snp)

getENSG <- function(rs = "rs180058", mart = mart.snp) {
  results <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"),
                   filters    = "snp_filter", values = rs, mart = mart)
  return(results)
}
genes <- map_dfc(snp_names, getENSG) %>% `[[`(ensembl_gene_stable_id)
edb <- query(AnnotationHub(),"EnsDb.Hsapiens.V90")
edb <- edb[[1]]
ensembldb::select(edb,keys = genes,columns = c("ENTREZID", "GENENAME"), keytype = "GENEID")
