fsdata <- readRDS("data\\fsdata_1104.rds")
dict_items <- dicGenerate(fsdata, exclude="rs|place|dat|sig|urotestd")
fsdata.new <- fsdata %>% mutate(stroke_subtype = str_replace_all(stroke_subtype,
                                                                 '[miss|other|unknown].+',
                                                                 'other'),
                                stroke_proven = ifelse(stroke_proven == "缺失", 
                                                        NA, 
                                                        stroke_proven),
                                ethnic2 = str_replace_all(ethnic2,
                                                         '族',
                                                         ''),
                                menaryear = ifelse(menoreason == "12",
                                                   12,
                                                   menaryear),
                                menoreason = ifelse(menoreason == "12",
                                                    "自然绝经",
                                                    menoreason),
                                ever_smoke = ifelse(smoke == "从不吸",
                                                    0,
                                                    1),
                                ever_drink = ifelse(drink == "从不饮酒",
                                                    0,
                                                    1),
                                druglipo = ifelse(drugaspr1 == "他汀类",
                                                  "是",
                                                  druglipo),
                                vitmin1 = ifelse(vitmin1 == "-6",
                                                  "否",
                                                 vitmin1),
                                qvsfs1 = ifelse(qvsfs1 == "9",
                                                "否",
                                                qvsfs1),
                                qvsfs2 = ifelse(qvsfs2 == "9",
                                                "否",
                                                qvsfs2),
                                qvsfs4 = ifelse(qvsfs4 == "9",
                                                "否",
                                                qvsfs4),
                                mostroke = ifelse(mostroke == 3,
                                                  "否",
                                                  mostroke),
                                mohemo = ifelse(mohemo == 0,
                                                  "否",
                                                mohemo),
                                facoro = ifelse(facoro == 9,
                                                "否",
                                                facoro),
                                mocoro = ifelse(mocoro == 3,
                                                "否",
                                                mocoro),
                                sibhemo = ifelse(sibhemo %in% c(0,5,9),
                                                "否",
                                                sibhemo),
                                hisstroke = ifelse(hisstroke == 0,
                                                   "否",
                                                   hisstroke),
                                memid = ifelse(memid == "父亲同胞1",
                                               "同胞1",
                                               memid))
write_rds(fsdata.new,"data\\fsdata_1106.rds")
write_excel_csv(fsdata.new,"data\\fsdata_1106_UTF8.csv")
## 
loadObj(fsdata.new, "data\\fsdata_1106.rds")
dict <- dicGenerate(fsdata.new)
write_excel_csv(dict,"audit\\dict.csv")

## 
dict <- read_delim("audit\\dict.txt", delim = "\t")
fsdata.label <- Value2Label(fsdata.new, dict)
write_rds(fsdata.label,"data\\fsdata_1112.rds")
write_excel_csv(fsdata.label,"data\\fsdata_1112_UTF8.csv")

## turn genotypes to sage formats

