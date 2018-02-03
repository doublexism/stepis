## hwe

hwe <- read_csv("audit/hwe.csv")
hwe_exclude <- read_csv("audit/HWE_excludeis.csv")
SNP_names <- str_subset(colnames(hwe),'^rs')

hwe_pvalue1 <- map_dfr(SNP_names, hwe_unrelate, data=hwe)
hwe_pvalue2 <- map_dfr(SNP_names, hwe_unrelate, data=hwe_exclude)

colnames(hwe_pvalue2) <- c("chi2", "pvalue2")
results <- bind_cols(hwe_pvalue1, hwe_pvalue2)
results$SNP <- SNP_names
results <- results %>% arrange(desc(chi)) %>% dplyr::select(SNP, chi, pvalue, chi2,pvalue2)
results[2:5] <- map(results[2:5], round, digits=4)
write_excel_csv(results, "audit/hwe.csv")
