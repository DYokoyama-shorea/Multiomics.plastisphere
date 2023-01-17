rm(list = ls())

library(tidyverse)
library(corrr)

CWD <- getwd()
DIR_Integrate <- str_c(CWD, "/res/integral") #"D:/Analysis/Analysis_Kikuchi_lab/2021_thesis_polymer.degradation"

Metabolite.2dj_output <- read_rds(str_c(DIR_Integrate, "/Metabolite.2dj_output.rds"))
#ProkAtlas_output <- read_rds(str_c(DIR_Integrate, "/ProkAtlas_output.rds"))
Microbiome_output <- read_rds(str_c(DIR_Integrate, "/Microbiome_output.rds"))
#PICRUSt2_output <- read_rds(str_c(DIR_Integrate, "/PICRUSt2_output.rds"))


Metabolite_2dj_org <- 
  Metabolite.2dj_output %>%
  filter(norm == "norm_pqn") %>%
  select(data) %>%
  unnest() %>%
  gather(-c(test, day, material, material2, material3), key = "mtb", value = "prop") %>%
  mutate(mtb = str_c("Metabolite_", mtb)) %>%
  spread(key = mtb, value = prop)

'''
Metabolite_2dj.annotated_org <-
  Metabolite_2dj_org %>%
  gather(-c(test, day, material, material2), key = "param", value = "prop") %>%
  mutate(param2 = str_sub(param, 12, -1)) %>% # Metabolie_の文字列を消す
  mutate(param3 = str_split(param2, pattern = "_")) %>% #_1,_2などの区分をなくすため
  mutate(num = map_int(param3, function(x){length(x)})) %>%　#分離した個数
  mutate(param4 = ifelse(!str_detect(param2, "ROI_") & num == 2, map(param3, function(x){x[1]}), param2)) %>%  #"ROI"がついていなくて、２つに分離(_1,_2とついていたもの)については、物質名を抽出
  unnest(param4) %>%
  group_by(test, day, material, material2, param4) %>%
  summarize(prop = sum(prop)) %>%
  filter(!str_detect(param4, "ROI")) %>%
  #ggplot(aes(x = param4, y = prop, color = material2)) +
  #geom_boxplot()
  spread(key = param4, value = prop)
'''

'''
NAME_2dj <-
  Metabolite_2dj_org %>%
  names() %>%
  as_tibble() %>%
  mutate(value = ifelse(str_detect(value, "Metabolite_"), str_sub(value, 12, -1), value)) %>%
  .$value
Metabolite_2dj.annotated_org <- #物質ごとにまとめない
  Metabolite_2dj_org %>% 
  set_names(NAME_2dj)
'''


Microbiome_org <- 
  Microbiome_output %>%
  filter(category == "g") %>%
  select(data4) %>%
  unnest()

Microbiome_top90_org <-
  Microbiome_org %>%
  mutate(id = row_number()) %>% # for adding id 
  gather(-c(test, material, material2, material3, day, bs, id), key = "taxa", value = "prop") %>%
  nest(-c(test, material, material2, day, bs, id)) %>%
  mutate(data2 = map(data, 
                     function(x){
                       x %>%
                         arrange(desc(prop)) %>%
                         mutate(cumsum = cumsum(prop)) %>%
                         filter(cumsum < 0.90)
                     })) %>%
  nest() %>%
  mutate(top90.list = map(data,
                          function(x){
                            x %>%
                              unnest(data2) %>%
                              .$taxa %>% 
                              unique()
                          })) %>%
  unnest(data) %>%
  mutate(data3 = map2(data, top90.list,
                      function(x, y){
                        x %>%
                          filter(taxa %in% y)
                      })) %>%
  select(-c(data, data2, top90.list)) %>%
  unnest(data3) %>%
  spread(key = taxa, value = prop) %>%
  select(-id)

Microbiome_over1p_org <-
  Microbiome_org %>%
  mutate(id = row_number()) %>% # for adding id 
  gather(-c(test, material, material2, material3, day, bs, id), key = "taxa", value = "prop") %>%
  nest(-c(taxa)) %>%
  mutate(max = map_dbl(data, function(x){max(x$prop)})) %>%
  filter(max > 0.01) %>%
  select(-max) %>%
  unnest(data) %>%
  spread(key = taxa, value = prop) %>%
  select(-id)

'''
PICRUSt2.ko_org <- 
  PICRUSt2_output %>%
  filter(db == "data_ko.2") %>%
  select(data3) %>%
  unnest() %>%
  select(-time)

PICRUSt2.path_org <- 
  PICRUSt2_output %>%
  filter(db == "data_pathway.2") %>%
  select(data3) %>%
  unnest() %>%
  select(-time) 


ProkAtlas_org <- 
  ProkAtlas_output %>% 
  select(category_2, data) %>% 
  unnest() %>% 
  spread(key = category_2, value = prop) %>% 
  select(-c(sample, time))
'''


####################################
# combine
####################################

DF_mba_mic <-
  Microbiome_over1p_org %>%
  left_join(., Metabolite_2dj_org, by = c("test", "day", "material", "material2", "material3")) %>%
  mutate(bs = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
 # na.omit() %>%
  select(-material) %>%
  mutate_if(is.character, as.factor) %>%
  select(-Unassigned) %>%
  nest(.key = "mic")
'''
DF_mba_ko <-
  PICRUSt2.ko_org %>%
  left_join(., Metabolite_2dj_org, by = c("test", "day", "material", "material2")) %>%
  na.omit() %>%
  select(-material) %>%
  mutate_if(is.character, as.factor) %>%
  nest(.key = "ko")
'''
DF_combine <-
  #bind_cols(DF_mba_mic, DF_mba_ko) %>%
  #gather(key = "param", value = "data")
  DF_mba_mic %>%
  gather(key = "param", value = "data")

write_rds(DF_mba_mic, str_c(DIR_Integrate, "/DF_mba_mic.rds"))
#write_rds(DF_mba_ko, str_c(DIR_Integrate, "/DF_mba_ko.rds"))
write_rds(DF_combine, str_c(DIR_Integrate, "/DF_combine.rds"))






