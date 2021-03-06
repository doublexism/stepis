loadObj(stk_dat, "data/stk_dat.rds")
loadObj(fsdata,"data/stroke_dat_0204.rds")
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
## discriptive family information
#table 1 family info
fsdata_summary <- setDT(fsdata)[,.(sib_num = sum(memid %in% c(1,4,5,6,7,8,9)), 
                            fa_num = sum(memid == 2), 
                            mo_num = sum(memid == 3)), family_no][,.(sib_num,fa_num,mo_num,pa_num = fa_num+mo_num)]

table1 <- fsdata_summary[,.(sum_pop = .N, sum_no_pa = sum(pa_num == 0),sum_fa = sum(fa_num), sum_mo = sum(mo_num), sum_pa = sum(pa_num ==2)), sib_num] %>%
  arrange(sib_num)

# table 2 disease status
fsdata_summary_dis <- setDT(fsdata)[memid %in% c(1,4,5,6,7,8,9),.(sib_num = .N, 
                                   case_only = (sum(stroke_proven) ==.N), 
                                   control_only = (sum(stroke_proven) == 0),
                                   discordant = (sum(stroke_proven) < .N & sum(stroke_proven)>0)) ,family_no]
table2 <- fsdata_summary_dis[,.(sum_co = sum(case_only), sum_con = sum(control_only), sum_dis = sum(discordant)), sib_num] %>%
  arrange(sib_num)

# discriptive analysis -- all
fsdata <- setDT(fsdata)
var_con <- c("age","bmi_calc","eduyear")
var_cat <- c("sex","edu","occu1","marriage","smoke2","drink2","dm_2",
             "ht","hyperlipidemia","hisatrial","hiscoronary","bmi_cat","cobesity")
con_desc <- map(fsdata[,var_con,with=FALSE],
                function(x, f) map_dbl(f, function(y) round(y(x, na.rm = TRUE),1)), 
                list(mean, sd)) %>% 
  do.call(rbind,.)
cat_desc <- map(fsdata[,var_cat,with=FALSE],
                table) 

con_desc_is <- map(fsdata[stroke_proven ==1,var_con,with=FALSE],
                   function(x, f) map_dbl(f, function(y) round(y(x, na.rm = TRUE),1)), 
                   list(mean, sd)) %>% 
  do.call(rbind,.)
con_desc_con <- map(fsdata[stroke_proven ==0,var_con,with=FALSE],
                   function(x, f) map_dbl(f, function(y) round(y(x, na.rm = TRUE),1)), 
                   list(mean, sd)) %>% 
  do.call(rbind,.)

cat_desc_is_con <- map(fsdata[,var_cat,with=FALSE],
                table,fsdata$stroke_proven == 1) 


## allele info
snp_names <- getSNPName(fsdata)
snp_names <- getSNPName(stk_dat)
## association test
geno_freq_table <- map_dfr(snp_names, genoFreq, data = stk_dat, Byvar = "stroke_proven")

# conditional logistic regression
snps <- colnames(stk_dat) %>% str_subset("^rs")
stk_add <- map(stk_dat[snp], )
cor_mat <- cor(stk_dat[snps])

library(survival)
clogit(is_stroke ~ rs2043211 + strata(fid), data = stk_dat,method = "exact") %>% summary()

