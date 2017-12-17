loadObj(stk_dat, "data/stroke_dat_1201.rds")
## association test

# conditional logistic regression
snps <- colnames(stk_dat) %>% str_subset("^rs")
stk_add <- map(stk_dat[snp], )
cor_mat <- cor(stk_dat[snps])

library(survival)
clogit(is_stroke ~ rs2043211 + strata(fid), data = stk_dat,method = "exact") %>% summary()

