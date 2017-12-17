select <- dplyr::select
#source("/Users/hsxmy/Desktop/fsdata/scripts/functions.R")
library(readr)
library(lubridate)
library(tidyverse)
fsdata <- read_rds("data/fsdata_1201.rds")

## diabetes
fsdata[c("pid","dm","hba1c2","fbg")][1:20,]

## hba1c2 >= 6.5 then diabetes
fsdata$dm <- ifelse(fsdata$hba1c2 >= 6.5, 1, fsdata$dm)

dm<-fsdata %>% filter(dm == 0, fbg > 7) %>% select(pid, dm, fbg, hba1c2, fasting)

## sex
fsdata[c("pid","sex","menopause","menoreason","preg")][1:20,]


## age missing data
fsdata$age <- ifelse(is.na(fsdata$age), 
                     interval(mdy(fsdata$birthday),
                     mdy(fsdata$survday))/years(1), fsdata$age)

fsdata$age <- ifelse(is.na(fsdata$age), 
              as.numeric(fsdata$survyear)- year(mdy(fsdata$birthday)), 
              fsdata$age)
agedf <- fsdata %>% filter(age <= 18)
pids_age <- agedf$pid
pids_age2 <- agedf$pid %>% setdiff("11040311")
fsdata <- fsdata %>% mutate(birthday = ifelse(pid %in% pids_age2, NA, birthday),
                            age = ifelse(pid %in% pids_age, NA, age) )

fsdata$age <- ifelse(is.na(fsdata$age), 
                     interval(mdy(fsdata$birthday),
                              mdy(fsdata$survday))/years(1), fsdata$age)
fsdata$age <- ifelse(is.na(fsdata$age), 
                     as.numeric(fsdata$survyear)- year(mdy(fsdata$birthday)), 
                     fsdata$age)
agedf <- fsdata %>% filter(is.na(fsdata$age))

## 
test <- fsdata %>% filter(pid %in% pids_age)

## occupation
meno <- fsdata %>% filter(!is.na(menopause)) %>% select(pid, menopause, menoage)
fsdata <- fsdata %>% mutate(menopause = ifelse(!is.na(menoage), 1, menopause))
fsdata[c("pid","menopause","menoage","menoreason")][1000:1020,]
meno <- fsdata %>% filter(sex == 1, !is.na(menopause) | !is.na(menoage) | !is.na(menoreason))
pids <- meno$pid
fsdata <- fsdata %>% mutate(sex = ifelse(pid %in% pids, 0, sex)) 

## smoke
fsdata <- fsdata %>% mutate(temp = ifelse(survyear %in% c(2013,2015) & smoke == 2,
                                          3,
                                          smoke),
                            temp = ifelse(survyear %in% c(2013,2015) & smoke == 3,
                                          2,
                                          temp))
                          
smoke <- fsdata %>% filter(temp == 1, !is.na(quitsmyear) | !is.na(smokeamout) | !is.na(smokeyear))
pid1 <- smoke$pid[1]
pid2 <- smoke$pid[2]
fsdata$temp[fsdata$pid == pid1] <- 2
fsdata$temp[fsdata$pid == pid2] <- 3
fsdata$smokeyear[fsdata$smokeyear == 0] <- NA
fsdata$quitsmyear[fsdata$quitsmyear == 0] <- NA

smoke <- fsdata %>% filter(temp == 3, !is.na(quitsmyear)) %>% select(pid, smoke, quitsmyear)
pid3 <- smoke$pid
fsdata <- fsdata %>% mutate(temp = ifelse(pid %in% pid3, 2, temp))
fsdata$smoke3 <- fsdata$temp

fsdata[["smoke2"]][fsdata$smoke3 == 1] <- 0
fsdata[["smoke2"]][fsdata$smoke3 %in% c(2,3)] <- 1
fsdata <- fsdata %>% select(-smoke, -temp)
table(fsdata$smoke3, !is.na(fsdata$quitsmyear), fsdata$survyear)
## drink
drink <- fsdata[c("pid","drink","drinkyear","quitdryear","drinkday","avedrink","maxdrink","onedrday")] %>% 
  filter(drink != 1) 
fsdata <- fsdata %>% mutate(temp = ifelse(survyear %in% c(2013,2015) & drink == 2,
                                         3,
                                         drink),
                           temp = ifelse(survyear %in% c(2013,2015) & drink == 3,
                                         2,
                                         temp),
                           drink3 = temp)
fsdata$quitdryear[fsdata$quitdryear == 0] <- NA
fsdata$drinkyear[fsdata$drinkyear == 0] <- NA
fsdata <- fsdata %>% mutate(drink3 = ifelse(!is.na(quitdryear), 2, drink3))
fsdata[["drink2"]][fsdata$drink3 == 1] <- 0
fsdata[["drink2"]][fsdata$drink3 %in% c(2,3)] <- 1
table(fsdata$drink2, fsdata$drink3)
fsdata <- fsdata %>% select(-temp, -drink)
fsdata.2015 <- filter(fsdata, survyear == "2015")

## stroke

fsdata$hishemo[fsdata$pid == "110610601"] <- 1
fsdata <- fsdata %>% mutate(temp = ifelse(survyear %in% c(2013) & hiscoronary == 1,
                                          0,
                                          hiscoronary),
                            temp = ifelse(survyear %in% c(2013) & hiscoronary == 0,
                                          1,
                                          temp),
                            temp = ifelse(!is.na(hiscorodate1), 1, temp),
                            hiscoronary = temp) %>% select(-temp)

## atrial
table(fsdata$fastroke, fsdata$survyear)

## height
test <- fsdata %>% filter(height < 130 | height > 200) %>% select(pid, height, height_bf, heightcm, weight, sex)
fsdata <- fsdata %>% mutate(height = ifelse(pid == "120405903", 148, height),
                            height = ifelse(pid == "121400404", 181, height),
                            height = ifelse(pid == "170910704", 162, height)) 

test <- fsdata %>% filter(is.na(height)) %>% select(pid, height, height_bf, heightcm, weight, sex)
fsdata <- fsdata %>% mutate(height = ifelse(is.na(height), heightcm, height))

## weight
test <- fsdata %>% filter(weight < 30 | weight > 150) %>% select(pid, weight, weight_bf, weightkg_x10, weight_from_bf,maxweight2)
weights <- test$weight_from_bf
weights[2] <- 75
fsdata$weight[fsdata$pid %in% test$pid] <- weights
fsdata$weight[fsdata$pid == "170911001"] <- 67
test <- fsdata %>% filter(is.na(weight)) %>% select(pid, weight, weight_bf, weight_from_bf, weightkg_x10, sex)

fsdata$weightkg_x10[fsdata$weightkg_x10 == 0] <- NA
fsdata <- fsdata %>% mutate(weight = ifelse(is.na(weight), weightkg_x10/10, weight))

## bmi
fsdata$bmi[fsdata$pid == "170410804"] <- NA
fsdata <- fsdata %>% mutate(bmi_calc = (weight/(height/100)**2))
fsdata$bmi_cat <- NA
fsdata$bmi_cat[fsdata$bmi_calc < 18.5] <- 1
fsdata$bmi_cat[fsdata$bmi_calc < 24 & fsdata$bmi_calc >= 18.5] <- 2
fsdata$bmi_cat[fsdata$bmi_calc < 28 & fsdata$bmi_calc >= 24] <- 3
fsdata$bmi_cat[fsdata$bmi_calc >= 28] <- 4


## bapwv and abi
fsdata <- fsdata %>% 
  rowwise() %>%
  mutate(abi = max(labi, rabi, na.rm = TRUE),
        bapwv = max(lbapwv, rbapwv, na.rm = TRUE))
fsdata$abi[fsdata$abi == -Inf] <- NA
fsdata$abi[fsdata$bapwv == Inf] <- NA
fsdata <- as_tibble(fsdata)
## CIMT
fsdata <- fsdata %>% mutate(cca_imt = (fimtavg_l1+fimtavg_l2+fimtavg_r1+fimtavg_r2)/4,
                            bif_imt = (fimtavg_r3 + fimtavg_l3)/2)


## drink raw
fs2015 <- read_csv("/Users/hsxmy/Desktop/fsdata/fs2015.csv")

table(fs2015$drink, !is.na(fs2015$quitdryear))
write_rds(fsdata, "/Users/hsxmy/Desktop/fsdata/stroke_dat_1201.rds")
### rsid
patch <- read_rds("/Users/hsxmy/Desktop/fsdata/fsdata.patch.rds")
fsdata <- fsdata %>% select(-rs556621,-rs2634074)

fsdata <- fsdata %>% left_join(patch, by = "pid")
write_rds(fsdata, "/Users/hsxmy/Desktop/fsdata/stroke_dat_1201.rds")
write_excel_csv(fsdata, "/Users/hsxmy/Desktop/fsdata/stroke_dat_1201.csv")
## maf
getAlleleFrequency(c("rs2634074"),data = fsdata)


##loadObj(stk, "data/stroke_dat_1201.rds")
write_excel_csv(fsdata,"data/fsdata_1215.csv")
