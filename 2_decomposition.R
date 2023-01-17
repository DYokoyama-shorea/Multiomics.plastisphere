rm(list = ls())

library(tidyverse)
library(lubridate)
library(scales)
library(readxl)
library(colorspace)
library(patchwork)

CWD <- getwd()
Dir_data_degradation <- str_c(CWD, "/data") # "D:/data/Miseq/Qiime2"
Dir_res_degradation <- str_c(CWD, "/res/degradation")
Dir_res_final_fig <- str_c(CWD, "/res/final.fig")

if(!file.exists(Dir_res_degradation)){
  dir.create(Dir_res_degradation)
}
if(!file.exists(Dir_res_final_fig)){
  dir.create(Dir_res_final_fig)
}


Order.material2 <-
  c("PHBH", "PCL", "PBSA", "PBS", "PBAT")

Order.material3 <-
  c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT")

df <- read_xlsx(str_c(Dir_data_degradation, "/degradation.xlsx"), sheet = 1)


df2 <-
  df %>%
  filter(wb == "water") %>%
  mutate(days = ifelse(str_detect(day, "w"), as.double(str_sub(day, 1,-2))*7,
                       ifelse(str_detect(day, "h"), as.double(str_sub(day, 1,-2))/24,
                              ifelse(str_detect(day, "d"), as.double(str_sub(day, 1,-2))*1, day)))) %>%
  type_convert() %>%
  mutate(material2 = ifelse(str_detect(material, "PCL"), "PCL",
                            ifelse(str_detect(material, "PHBH"), "PHBH",
                                   ifelse(str_detect(material, "PBSA"), "PBSA",
                                          ifelse(str_detect(material, "PBS"), "PBS",
                                                 ifelse(str_detect(material, "PBAT"), "PBAT", material)))))) %>%
  mutate(material2 = factor(material2, levels = Order.material2)) %>%
  mutate(material3 = ifelse(str_detect(material, "PCL"), "PCL",
                            ifelse(str_detect(material, "PHBH6"), "PHBH_HH6%",
                                   ifelse(str_detect(material, "PHBH"), "PHBH_HH10%",
                                          ifelse(str_detect(material, "PBSA"), "PBSA",
                                                 ifelse(str_detect(material, "PBS"), "PBS",
                                                        ifelse(str_detect(material, "PBAT"), "PBAT", material))))))) %>%
  mutate(material3 = factor(material3, levels = Order.material3)) %>%
  mutate(g = str_c(test, material, sep = "_")) %>%
  mutate(test = ifelse(test == "test6", "Test.A",
                       ifelse(test == "test8", "Test.B",
                              ifelse(test == "test9", "Test.C", test)))) 


DF_decompose <-
  df2 %>%
  mutate(ratio.decompose = 100 * ratio.decompose) %>%
  gather(decrease, ratio.decompose, `BOD-ctrl`, key = "method", value = "value") %>%
  nest(-method) %>%
  mutate(ylab = ifelse(method == "decrease", "Weight Loss (mg)",
                       ifelse(method == "ratio.decompose", "Weight Loss (%)",
                              ifelse(method == "BOD-ctrl", "BOD", NA)))) %>%
  mutate(plot = map2(data, ylab,
                    function(x,y){
                      x %>%
                        ggplot(aes(x = days, y = value, color = material3, group = g)) +
                        geom_point() +
                        geom_line() +
                        facet_grid(test ~ .)+
                        theme_bw() +
                        labs(y = y, x = "Incubation Day", color = "Material") +
                        scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                      "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                      "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                      "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                      "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                      "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
                        
                    }))

df_bod <- read_xlsx(str_c(Dir_data_degradation, "/degradation.xlsx"), sheet = "BOD")

plot_bod <-
  df_bod %>%
  gather(-day, key = "key", value = "value") %>%
  separate(key, into = c("test", "material"), sep = "_") %>%
  filter(material != "control") %>%
  mutate(material2 = ifelse(str_detect(material, "PCL"), "PCL",
                            ifelse(str_detect(material, "PHBH"), "PHBH",
                                   ifelse(str_detect(material, "PBSA"), "PBSA",
                                          ifelse(str_detect(material, "PBS"), "PBS",
                                                 ifelse(str_detect(material, "PBAT"), "PBAT", material)))))) %>%
  mutate(material2 = factor(material2, levels = Order.material2)) %>%
  mutate(material3 = ifelse(str_detect(material, "PCL"), "PCL",
                            ifelse(str_detect(material, "PHBH6"), "PHBH_HH6%",
                                   ifelse(str_detect(material, "PHBH"), "PHBH_HH10%",
                                          ifelse(str_detect(material, "PBSA"), "PBSA",
                                                 ifelse(str_detect(material, "PBS"), "PBS",
                                                        ifelse(str_detect(material, "PBAT"), "PBAT", material))))))) %>%
  mutate(material3 = factor(material3, levels = Order.material3)) %>%
  na.omit() %>%
  ggplot(aes(x = day, y = value, color = material3)) +
  #geom_point(alpha = 0.6) + 
  geom_line(alpha = 0.6, size = 1)  +
  facet_grid(test ~ .) + 
  theme_bw() +
  labs(x = "Incubation Day", y = "BOD consumption (mg/L)", color = "Material") +
  scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))


ggsave(plot = plot_bod, file = str_c(Dir_res_degradation, "/BOD_consumption.png", sep = ""), height = 4, width = 4.5)
for(i in 1:(nrow(DF_decompose))){
  plot_tmp <- DF_decompose[i,4] %>% .[[1]] %>% .[[1]]
  ggsave(plot = plot_tmp, file = str_c(Dir_res_degradation, "/", DF_decompose[i,1] %>% .[[1]], ".png", sep = ""), height = 4, width = 4.5)
}
  


df_dna <- 
  read_xlsx(str_c(Dir_data_degradation, "/DNA.xlsx"), sheet = 1) %>%
  mutate(DNA_surf_ng.ml = ifelse(str_detect(DNA_surf_ng.ml, "<|-"), NA, DNA_surf_ng.ml),
         DNA_bio_ng.ml = ifelse(str_detect(DNA_bio_ng.ml, "<|-"), NA, DNA_bio_ng.ml)) %>%
  mutate(days = ifelse(str_detect(day, "w"), as.double(str_sub(day, 1,-2))*7,
                       ifelse(str_detect(day, "h"), as.double(str_sub(day, 1,-2))/24,
                              ifelse(str_detect(day, "d"), as.double(str_sub(day, 1,-2))*1, day)))) %>%
  type_convert() %>%
  mutate(material2 = ifelse(str_detect(material, "PCL"), "PCL",
                            ifelse(str_detect(material, "PHBH"), "PHBH",
                                   ifelse(str_detect(material, "PBSA"), "PBSA",
                                          ifelse(str_detect(material, "PBS"), "PBS",
                                                 ifelse(str_detect(material, "PBAT"), "PBAT", material)))))) %>%
  mutate(material2 = factor(material2, levels = Order.material2)) %>%
  mutate(material3 = ifelse(str_detect(material, "PCL"), "PCL",
                            ifelse(str_detect(material, "PHBH6"), "PHBH_HH6%",
                                   ifelse(str_detect(material, "PHBH"), "PHBH_HH10%",
                                          ifelse(str_detect(material, "PBSA"), "PBSA",
                                                 ifelse(str_detect(material, "PBS"), "PBS",
                                                        ifelse(str_detect(material, "PBAT"), "PBAT", material))))))) %>%
  mutate(material3 = factor(material3, levels = Order.material3)) %>%
  mutate(test = ifelse(test == "test6", "Test.A",
                       ifelse(test == "test8", "Test.B",
                              ifelse(test == "test9", "Test.C", 
                                     ifelse(test == "test10", "Test.D", 
                                            ifelse(test == "test12", "Test.E", test)))))) 

plot.dna <-
  df_dna %>%
  filter(test == "Test.D" | test == "Test.E") %>%
  mutate(DNA_surf_ng.cm2 = DNA_surf_ng.ml * vol_surf_ml / area_cm2, 
         DNA_bio_ng.cm2 = DNA_bio_ng.ml * vol_bio_ml / area_cm2) %>%
  #mutate(DNA_ng.cm2 = DNA_surf_ng.cm2 + DNA_bio_ng.cm2) %>%
  mutate(id = str_c(test, material)) %>%
  gather(DNA_surf_ng.cm2, DNA_bio_ng.cm2, key = "bs", value = "DNA_ng.cm2") %>%
  mutate(bs = ifelse(bs == "DNA_surf_ng.cm2", "Surface", "Biofilm")) %>%
  #group_by(test, days, material, material2, id) %>%
  #summarize(DNA_ng.cm2 = sum(DNA_ng.cm2, na.rm = T)) %>%
  mutate(DNA_ug.cm2 = DNA_ng.cm2/1000) %>%
  ggplot(aes(x = days, y = DNA_ug.cm2, color = material3, group = id)) +
  geom_point() +
  geom_line() +
  facet_grid(test ~ bs, scales = "free") +
  theme_bw() +
  labs(x = "Incubation Days", y = "DNA (ug/cm2)", color = "Material") +
  scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7]#, 
                                #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]
                                )) 



######################
### Figs for thesis
######################

Figure.S2 <-  DF_decompose[[2,4]][[1]] + theme(legend.position = "none") + plot_bod
ggsave(plot = Figure.S2, file = str_c(Dir_res_final_fig, "/Figure.S2.pdf", sep = ""), height = 4, width = 8.5)

#ggsave(plot = plot.dna, file = str_c(Dir_res_final_fig, "/Figure.S3.pdf", sep = ""), height = 4, width = 8)
