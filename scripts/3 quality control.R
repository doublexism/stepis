stk_dat <- read_rds("data\\stk_dat.rds")
## hardy weinburg test
# Based on family data
# 1. sample siblings from each family
sib_code <- c(1,4:9)
hwe_fam <- stk_dat %>% 
  group_by(family_no) %>% 
  mutate(n_sib = sum(memid %in% sib_code)) %>%
  filter(n_sib >= 2) %>%
  mutate(ind = rank(family_no, ties.method = "last")) %>%
  filter(ind %in% c(1,2))
# 2. clean data
SNP_names <- getSNPName(hwe_fam)
hwe_fam[SNP_names] <- map(hwe_fam[SNP_names], geno_mat)
# 2. hwe test
hwe_sib <- map_dfr(SNP_names, hwe_test, data = hwe_fam, famid = "family_no")
hwe_sib[2:3] <- map(hwe_sib[2:3], round, digits = 4) 
hwe_sib <- hwe_sib %>% arrange(p)
# 3. write to csv
write_excel_csv(hwe_sib, "data\\hwe.csv")

## output to sage formate
SNP_names <- getSNPName(stk_dat)
vars_sage <- c("family_no","pid","fatherid","motherid","sex","stroke_proven","age", "ever_smoke", 
               "ever_drink", "dm", "ht", SNP_names)
sage_dat <- stk_dat %>% 
  select(one_of(vars_sage)) %>%
  mutate_at(vars(one_of(SNP_names)), sageSNP) %>%
  mutate_all(sageVar)
# genomic info files
# SNP_info <- Get.SNPinfo(SNP_names, "CHB")
# write_rds(SNP_info, "sage\\snpinfo.rds")
loadObj(SNP_info, "sage\\snpinfo.rds")
SNP_info[c("chromesome","site")] <- map(SNP_info[c("chromesome","site")], as.numeric)
sage_gmf <- SNP_info %>% 
  mutate(rsid = paste0(chromesome,"_",SNP)) %>%
  arrange(chromesome, site) %>%
  select(rsid, site)
# output to sage
write_csv(sage_dat, "sage\\sage_1114.csv")
write_csv(sage_gmf, "sage\\gmf_1116.csv")
## qc 
error_data <- stk_dat %>% group_by(family_no) %>% 
  mutate(np = sum(memid %in% c(2,3)),
         page = min(ifelse(memid %in% c(2,3), age, NA), na.rm = TRUE)) %>%
  select(pid, family_no, name,memid,sex, age,birthday,id, menoage,menopause,menoreason,sex_abi, sex_bf, page, np)%>%
  #filter(np >= 1 & (page - age < 15) & memid == 1 )
  filter(family_no %in% c("1421107", "1709117", "1096","1702107"))

