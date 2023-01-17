rm(list = ls())


###
#2dj
###

library(readxl)
library(tidyverse)
library(Rcpm)
library(vegan)
library(ggrepel)
library(pairwiseAdonis)
library(colorspace)
library(ggpubr)
library(patchwork)

CWD <- getwd()
Dir_data <- str_c(CWD, "/data", sep = "")
Dir_res_2dj <- str_c(CWD, "/res/2dj")
Dir_res_integral <- str_c(CWD, "/res/integral")
Dir_res_final_fig <- str_c(CWD, "/res/final.fig")

if(!file.exists(Dir_res_2dj)){
  dir.create(Dir_res_2dj)
}
if(!file.exists(Dir_res_integral)){
  dir.create(Dir_res_integral)
}
if(!file.exists(Dir_res_final_fig)){
  dir.create(Dir_res_final_fig)
}

Order.material2 <-
  c("PHBH", "PCL", "PBSA", "PBS", "PBAT")
Order.material3 <-
  c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT")

Surface.area <-
  tibble(test = c("Test.A", "Test.B", "Test.C"),
         surface.area = c((1*6)*2 + ((1+6)*2)*0.05, (1*6)*2 + ((1+6)*2)*0.02, (0.5*6)*2 + ((0.5+6)*2)*0.02))

DF_2dj_raw_org <- 
  read_xlsx(str_c(Dir_data, "/2dj_roiSummary_t6-9.V2.xlsx", sep =""), skip = 1, sheet = 1) %>%
  mutate(time_num = str_sub(time, 1, -2) %>% as.double()) %>%
  mutate(day = ifelse(str_detect(time, "h"), time_num /24, 
                      ifelse(str_detect(time, "d"), time_num,
                             ifelse(str_detect(time, "w"), time_num*7, NA)))) %>%
  mutate(material2 = ifelse(str_detect(material, "PBSA"), "PBSA",
                            ifelse(str_detect(material, "PCL"), "PCL",
                                   ifelse(str_detect(material, "PHBH"), "PHBH", material)))) %>%
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
                                            ifelse(test == "test12", "Test.E", "Other"
                                            )))))) %>%
  dplyr::select(-c(num, time, time_num, File, DSS)) %>%
  dplyr::select(test, day, material, material2, material3, point, everything()) %>%
  filter(point == "water") %>% #指定しなければ(コメントアウトすれば)、waterとbiofilmを足した値を返す。
  dplyr::select(-point) %>%
  group_by(test, day, material, material2, material3)  %>%
  summarize_all(sum) %>%
  ungroup()

DF_2dj_raw.grouped_org <- #物質ごとにグループ化したデータ
  DF_2dj_raw_org %>%
  gather(-c(test, day, material, material2, material3), key = "ROI", value = "intensity") %>%
  mutate(ROI.2 = str_split(ROI, pattern = "_")) %>%
  mutate(num = map_int(ROI.2, function(x){length(x)})) %>%
  mutate(ROI.3 = ifelse(!str_detect(ROI, "ROI_") & num == 2, map(ROI.2, function(x){x[1]}), ROI)) %>%
  unnest(ROI.3) %>%
  group_by(test, day, material, material2, material3, ROI.3) %>%
  summarize(intensity = sum(intensity)) %>%
  spread(key = ROI.3, value = intensity) %>%
  ungroup()
  


DF_2dj <-
  nest(DF_2dj_raw_org, .key = "raw") %>% #生のROIを使いたければ、DF_2dj_raw_org.グループ分けしたデータを使いたければDF_2dj_raw.grouped_org
  mutate(raw_per_surf = map(raw,
                            function(x){
                              x %>%
                                left_join(., Surface.area, by = "test") %>%
                                mutate_at(vars(-c(test, day, material, material2, material3)), ~./surface.area) %>%
                                dplyr::select(-surface.area)
                            })) %>%
  mutate(norm_area = map(raw, 
                    function(x){
                      tmp <- x %>% dplyr::select(-c(test, day, material, material2, material3)) 
                      x %>% 
                        bind_cols(., tibble(sum = rowSums(tmp))) %>%
                        mutate_at(vars(-c(test, day, material, material2, material3)), ~./sum) %>%
                        dplyr::select(-sum)})) %>%
  mutate(norm_pqn  = map(norm_area, 
                    function(x){
                      tmp <- x %>% dplyr::select(-c(test, day, material, material2, material3)) %>% pqn() %>% as_tibble()
                      x %>% 
                        dplyr::select(c(test, day, material, material2, material3)) %>%
                        bind_cols(., tmp)})) %>%
  gather(key = "norm", value = "data") %>%
  mutate(nmds = map(data,
                    function(x){
                      x %>%
                        dplyr::select(-c(test, day, material, material2, material3)) %>%
                        metaMDS(k = 2, trymax = 20)
                    })) %>%
  mutate(plot.nmds = map2(data, nmds,
                         function(x, y){
                           x %>%
                             dplyr::select(c(test, day, material, material2, material3)) %>%
                             bind_cols(., y$points %>% as_tibble()) %>%
                             ggplot(aes(x = MDS1, y =MDS2, color = material3, fill = material3, size = day, shape = test)) +
                             geom_point(alpha = 0.5) +
                             scale_shape_manual(values = c("Test.A" = 21, "Test.B" = 22, "Test.C" = 23)) + # "test.D" = 24, "test.E" = 25
                             theme_bw()  +
                             theme(legend.direction = "vertical", legend.box = "horizontal") +
                             labs(title = str_c("Stress = ", round(y$stress, 3)),
                                  size = "Day", color = "Material", shape = "Test") +
                             guides(shape = guide_legend(order = 1),
                                    color = guide_legend(order = 2),
                                    size = guide_legend(order = 3),
                                    fill = F) +
                             scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                           "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                           "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                           "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                           "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                           "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                             scale_fill_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                          "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                          "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                          "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                          "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                          "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))
                         })) %>%
  mutate(plot.score = map(nmds,
                          function(x){
                            x$species %>% 
                              as_tibble(rownames = "ROI") %>%
                              mutate(annotated = ifelse(str_detect(ROI, "ROI"), "unannotated", "annotated"),
                                     label = ifelse(str_detect(ROI, "ROI"), NA_character_, ROI)) %>%
                              ggplot(aes(x = MDS1, y = MDS2, label = label, color = annotated)) +
                              geom_point() +
                              geom_text_repel(color = "#F8766D", size = 3) +
                              scale_color_manual(values = c("unannotated" = "grey", "annotated" = "black")) +
                              theme_bw() #+
                              #theme(legend.position = "none") 
                          })) %>%
  mutate(permanova = map(data,
                         function(x){
                           tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day)) 
                           tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day))
                           res <- adonis2(tmp1 ~ test*material3*day, data = tmp2, method = "bray")
                           return(res)
                         })) %>%
  mutate(pairwise.permanova = map(data,
                                  function(x){
                                    tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day)) 
                                    tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day))
                                    res <- pairwise.adonis(tmp1,tmp2$material3)
                                    return(res)
                                  })) %>%
#  mutate(pairwise.permanova.2 = map(data,
#                                  function(x){
#                                    tmp1 <- x %>% dplyr::select(-c(test, material, material2, day))  %>% as.data.frame()
#                                    tmp2 <- x %>% dplyr::select(c(test, material, material2, day)) %>% as.data.frame()
#                                    res <- pairwise.adonis(tmp1 ~ material2 * test, data = tmp2, strata = "material2")
#                                    return(res)
#                                  }))
  mutate(boxplot.all = map(data,
                           function(x){
                             x %>%
                               gather(-c(test, day, material, material2, material3), key = "ROI", value = "intensity") %>%
                               nest(-ROI) %>%
                               mutate(boxplot.each = map2(data, ROI,
                                                         function(x, y){
                                                           x %>%
                                                             ggplot(aes(x = material2, y = intensity, fill = material3, color = material3)) +
                                                             geom_boxplot(alpha = 0.6) +
                                                             theme_bw() +
                                                             scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                                           "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                                           "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                                           "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                                           "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                                           "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                                             scale_fill_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                                          "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                                          "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                                          "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                                          "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                                          "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))+
                                                             xlab("Material") +
                                                             ylab("Intensity") +
                                                             facet_grid(test ~ .) +
                                                             ggtitle(y)})) %>%
                               mutate(smooth.each = map2(data, ROI,
                                                         function(x, y){
                                                           x %>%
                                                             mutate(id = str_c(test, material, sep = "_")) %>%
                                                             ggplot(aes(x = day, y = intensity, group = id, color = material3)) +
                                                             geom_point() +
                                                             geom_line() +
                                                             #geom_smooth(se = F) +
                                                             theme_bw() +
                                                             scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                                           "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                                           "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                                           "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                                           "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                                           "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                                             labs(y = "Intensity", x = "Incubation Day", color = "Material") +
                                                             ggtitle(y)  
                                                         }))
                           })) 

for(i in 1:nrow(DF_2dj)){
  norm <- DF_2dj[i,] %>% dplyr::select(norm) %>% .[[1]]
  plot_nmds <- DF_2dj[i,] %>% dplyr::select(plot.nmds) %>% .[[1]] %>% .[[1]]
  plot_score <- DF_2dj[i,] %>% dplyr::select(plot.score) %>% .[[1]] %>% .[[1]]
  ggsave(plot = plot_nmds, file = str_c(Dir_res_2dj, "/nMDS_2dj_", norm, ".png", sep = ""), height = 3.5, width = 6.8)
  ggsave(plot = plot_score, file = str_c(Dir_res_2dj, "/nMDS.score_2dj_", norm, ".png", sep = ""), height = 3, width = 4.5)
}

write_rds(DF_2dj, str_c(Dir_res_integral, "/Metabolite.2dj_output.rds"))

#Final figure
Figure.S11_1 <- DF_2dj[[4,4]][[1]] + labs(x = "nMDS1", y = "nMDS2") +theme(legend.position = "none")
legend <- DF_2dj[[4,4]][[1]] %>% get_legend() %>% as_ggplot()
Figure.S11_2.1 <- DF_2dj[[4,8]][[1]][[1,3]][[1]] + theme(legend.position = "none")
Figure.S11_2.2 <- DF_2dj[[4,8]][[1]][[2,3]][[1]] + theme(legend.position = "none")
Figure.S11_2.3 <- DF_2dj[[4,8]][[1]][[6,3]][[1]] + theme(legend.position = "none")

Figure.S11 <- (Figure.S11_1 | legend)/(Figure.S11_2.1 | Figure.S11_2.2 | Figure.S11_2.3) + plot_layout(ncol = 1)
ggsave(plot = Figure.S11, file = str_c(Dir_res_2dj, "/Figure.2dj.pdf", sep = ""), height = 8, width = 8)


######################
### Figs for thesis
######################
ggsave(plot = Figure.S11, file = str_c(Dir_res_final_fig, "/Figure.S11.pdf", sep = ""), height = 8, width = 8)

