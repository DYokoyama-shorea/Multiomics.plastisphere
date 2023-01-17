rm(list = ls())

library(tidyverse)
library(lubridate)
library(scales)
library(vegan)
library(pairwiseAdonis)
library(ggh4x) # for extended facet in ggplot2
library(plotly)
library(ggpubr)
library(patchwork)
library(colorspace)
library(tictoc)
library(multcomp)
library(nparcomp)

library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("helper_functions.r")

CWD <- getwd()
Dir_data_miseq <- str_c(CWD, "/data") # "D:/data/Miseq/Qiime2"
Dir_res_microbiome <- str_c(CWD, "/res/microbiome")
Dir_res_microbiome_rarefaction <- str_c(Dir_res_microbiome, "/rarefaction")
Dir_res_integral <- str_c(CWD, "/res/integral")
Dir_res_final_fig <- str_c(CWD, "/res/final.fig")

if(!file.exists(Dir_res_microbiome)){
  dir.create(Dir_res_microbiome)
}
if(!file.exists(Dir_res_microbiome_rarefaction)){
  dir.create(Dir_res_microbiome_rarefaction)
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

DF3 <- 
  read_rds(str_c(Dir_data_miseq, "/microbe.rds", sep = ""))

DF3_rarefied <-
  read_rds(str_c(Dir_data_miseq, "/microbe_rarefied.rds", sep = "")) %>%
  mutate(data.rarefied.2 = map2(info, data.rarefied,
                                function(x, y){
                                  y %>%
                                    left_join(., x %>% rename(sample = `#SampleID`), by = "sample") %>%
                                    separate(Description, into = c("test", "description"), sep = ";") %>%
                                    separate(description, into = c("material", "time", "bs")) %>%
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
                                    ungroup() %>%
                                    dplyr::select(-c(sample, Project, time, time_num)) %>%
                                    filter(material2 != "water") %>%
                                    dplyr::select(test, material, material2, day, bs, everything()) %>%
                                    filter(material != "water") %>%
                                    filter(!(test == "test4" | test == "test5" | test == "test5S")) %>%
                                    mutate(test = ifelse(test == "test6", "Test.A",
                                                         ifelse(test == "test8", "Test.B",
                                                                ifelse(test == "test9", "Test.C",
                                                                       ifelse(test == "test10", "Test.D",
                                                                              ifelse(test == "test12", "Test.E", "Other"
                                                                              )))))) %>%
                                    dplyr::select(test, material, material2, material3, day, bs, everything())
                                })) %>%
  mutate(res_nmds = map(data.rarefied.2,
                        function(x){
                          set.seed(1)
                          tmp <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs))
                          res <- metaMDS(tmp, k = 2, trymax = 20)
                          return(res)
                        })) %>%
  mutate(plot_nmds = map2(data.rarefied.2, res_nmds,
                          function(x, y){
                            tmp <- x %>% dplyr::select(c(test, material, material2, material3, day, bs))
                            tmp2 <- y$points %>% as_tibble()
                            bind_cols(tmp, tmp2) %>%
                              ggplot(aes(x = MDS1, y = MDS2, color = material3, fill = material3, size = day, shape = test)) +
                              geom_point(alpha = 0.5) +
                              theme_bw() +
                              theme(legend.direction = "vertical", legend.box = "horizontal") +
                              labs(title = str_c("Stress = ", round(y$stress, 3)),
                                   size = "Day", color = "Material", shape = "Test") +
                              guides(shape = guide_legend(order = 1),
                                     color = guide_legend(order = 2),
                                     size = guide_legend(order = 3),
                                     fill = "none") +
                              scale_shape_manual(values = c("Test.A" = 21, "Test.B" = 22, "Test.C" = 23, "Test.D" = 24, "Test.E" = 25)) +
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
  mutate(permanova = map(data.rarefied.2,
                         function(x){
                           tic()
                           tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
                           tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
                           res <- adonis2(tmp1 ~ test + material3 +  day + bs, data = tmp2, permutations = 999, method = "bray")
                           toc()
                         return(res)
                         })) %>%
  mutate(permanova.2 = map(data.rarefied.2,
                         function(x){
                           tic()
                           tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
                           tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
                           res <- adonis2(tmp1 ~ test * material3 * day * bs, data = tmp2, permutations = 999, method = "bray")
                           toc()
                           return(res)
                         })) %>%
  mutate(data.bray.distance = map(data.rarefied.2,
                                  function(x){
                                    tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs))
                                    tmp2 <- x %>% dplyr::select(c(test, material3, day, bs)) %>% mutate(id = row_number() %>% as.character())
                                    tmp3 <- tmp1 %>% vegdist(method = "bray")
                                    tmp4 <-
                                      tmp3 %>%
                                      as.matrix() %>% as_tibble() %>%
                                      bind_cols(tmp2, .) %>%
                                      gather(-c(test,  material3, day, bs, id), key = "id2", value = "distance") %>%
                                      left_join(., tmp2 %>% set_names(c("test.to", "material.to", "day.to", "bs.to", "id2"))) %>%
                                      filter(test == test.to) %>%
                                      filter(id != id2) #%>%
                                    #mutate(pair = map2(material3, material.to, function(x,y){c(as.character(x),as.character(y))%>% sort()})) %>%
                                    #mutate(pair.1 = map_chr(pair, ~.[1]), pair.2 = map_chr(pair, ~.[2])) %>%
                                    #mutate(pair = str_c(pair.1, pair.2, sep = "_vs_"))
                                    
                                  })) %>%
  mutate(plot.bray.distance = map(data.bray.distance,
                                  function(x){
                                    tmp1 <-
                                      x %>%
                                      nest(-material.to, -test) %>%
                                      mutate(multcomp = map(data,
                                                            function(x){
                                                              #tic()
                                                              res <- nparcomp(distance ~ material3, data=x, asy.method = "mult.t",
                                                                              type = "Tukey",alternative = "two.sided",info = FALSE)
                                                              #toc()
                                                              return(res)
                                                            })) %>%
                                      mutate(letters = map(multcomp, func_letter_nparcomp))
                                    
                                    data.letters <- 
                                      tmp1 %>%
                                      dplyr::select(test, material.to, letters) %>%
                                      unnest() %>%
                                      rename(material3 = "param") %>%
                                      mutate(material3 = factor(material3, levels = Order.material3))
                                    
                                    data.distance <-
                                      tmp1 %>%
                                      dplyr::select(test, material.to, data) %>%
                                      unnest() %>%
                                      mutate(material3 = factor(material3, levels = Order.material3))
                                    
                                    ggplot() +
                                      geom_jitter(data = data.distance, aes(x = material3, y = distance, color = material3), alpha = 0.05) +
                                      geom_boxplot(data = data.distance, aes(x = material3, y = distance), fill = NA, outlier.colour = NA) +
                                      geom_text(data = data.letters, aes(x = material3, y = 1.15, label = fin.letter)) +
                                      theme_bw() +  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
                                      facet_grid(test ~ material.to) +
                                      ylim(0, 1.2) +
                                      ylab("Bray-Curtis Index") + xlab("Material") +
                                      scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                    "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                    "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                    "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                    "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                    "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                      theme(legend.position = "none")
                                    
                                    
                                  }))

for(i in 1:nrow(DF3_rarefied)){
  plot1 <- DF3_rarefied[i,] %>% dplyr::select(plot_nmds) %>% .[[1]] %>% .[[1]]
  plot2 <- DF3_rarefied[i,] %>% dplyr::select(plot.bray.distance) %>% .[[1]] %>% .[[1]]
  ggsave(plot = plot1, file = str_c(Dir_res_microbiome_rarefaction, "/nMDS_rarefied.", i, ".png"), width = 7, height = 4)
  ggsave(plot = plot2, file = str_c(Dir_res_microbiome_rarefaction, "/bray.distance_rarefied", i, ".png"), width = 12, height = 6)

}

write_rds(DF3_rarefied %>% dplyr::select(data.rarefied.2, res_nmds, permanova, permanova.2),
          str_c(Dir_res_microbiome, "/df_rarefied.rds"))
  
# 
# 
# res <- list()
# for(i in 1:nrow(DF3_rarefied)){
#   tic()
#   tmp <- DF3_rarefied[i,] %>% dplyr::select(data.rarefied.2) %>% unnest()
#   tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#   tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
#   tmp3 <- pairwise.adonis2(tmp1 ~ material3 + test + day + bs , data = tmp2)
#   res <- append(res, list(tmp3))
#   toc()
# }
# 
# res2 <- list()
# for(i in 1:nrow(DF3_rarefied)){
#   tic()
#   tmp <- DF3_rarefied[i,] %>% dplyr::select(data.rarefied.2) %>% unnest()
#   tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#   tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
#   tmp3 <- pairwise.adonis2(tmp1 ~ material3 * test * day * bs , data = tmp2)
#   res2 <- append(res2, list(tmp3))
#   toc()
# }
# 
# 
# DF3_rarefied_2 <-
#   bind_cols(DF3_rarefied, tibble(pairwise.permanova = res))  %>%
#   bind_cols(., tibble(pairwise.permanova.2 = res2))  %>%
#   mutate(permanova = map(permanova, function(x)(as_tibble(x, rownames = "factor"))))  %>%
#   mutate(permanova.2 = map(permanova.2, function(x)(as_tibble(x, rownames = "factor"))))  %>%
#   mutate(pairwise.permanova = map(pairwise.permanova, function(x){tibble(pair = names(x), res = x)})) %>%
#   mutate(pairwise.permanova.2 = map(pairwise.permanova.2, function(x){tibble(pair = names(x), res = x)})) 
# 
# save.tmp <-
#   DF3_rarefied_2 %>%
#   dplyr::select(#res_nmds, plot_nmds, 
#          permanova, permanova.2, pairwise.permanova, pairwise.permanova.2) 
# write_rds(save.tmp, "rarefaction_permanova.rds")
# 
# DF3_rarefied_2 %>% dplyr::select(permanova.2) %>% mutate(rarefaction = row_number()) %>% unnest() %>%
#   filter(factor != "Total") %>%
#   nest(-factor) %>%
#   mutate(factor = str_replace(factor, "3", "")) %>%
#   mutate(factor = str_replace(factor, ":", " x ")) %>%
#   mutate(factor = str_replace(factor, ":", " x ")) %>%
#   mutate(factor = str_replace(factor, ":", " x ")) %>%
#   mutate(factor= factor(factor, levels = .$factor)) %>%
#   unnest() %>%
#   ggplot(aes(x = factor, y = R2)) +
#   geom_boxplot() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# res3 <- tibble()
# for(i in 1:nrow(DF3_rarefied)){
#   TMP <- DF3_rarefied[i,] %>% dplyr::select(data.rarefied.2) %>% unnest() %>% nest(-test) %>% arrange(test)
#   for(j in 1:nrow(TMP)){
#     tic()
#     tmp <- TMP[j,] %>% unnest()
#     tmp.test <- TMP[j,] %>% dplyr::select(test) %>% .[[1]]
#     tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs))
#     tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs))
#     tmp3 <- try(pairwise.adonis2(tmp1 ~ material3 * day * bs , data = tmp2))
#     tmp4 <- try(pairwise.adonis2(tmp1 ~ material3 * day, data = tmp2))
#     tmp3 <- ifelse(tmp3 == "try-error", "NA", tmp3)
#     tmp4 <- ifelse(tmp4 == "try-error", "NA", tmp4)
#     res.tmp <- tibble(rarefaction = i, test = tmp.test, pairwise.permanova.1 = list(tmp3), pairwise.permanova.2 = list(tmp4))
#     res3 <- bind_rows(res3, res.tmp)
#     rm(tmp, tmp.test, tmp1, tmp2, tmp3, tmp4, res.tmp)
#     toc()
#   }
# }
# res3 %>%
#   dplyr::select(-pairwise.permanova.2) %>%
#   mutate(class = map_chr(pairwise.permanova.1, function(x){class(x) %>% .[[1]]})) %>%
#   filter(class != "character") %>%
#   mutate(data = map(pairwise.permanova.1, function(x){tibble(pair = names(x), res = x)})) %>%
#   dplyr::select(rarefaction, test, data) %>%
#   unnest() %>%
#   filter(pair != "parent_call") %>%
#   mutate(res = map(res, function(x){as_tibble(x, rownames = "factor")})) %>%
#   unnest() %>%
#   filter(factor == "material3") %>%
#   group_by(pair) %>%
#   summarize(R2_mean = mean(R2), R2_sd = sd(R2)) %>%
#   mutate(pair2 = str_split(pair, pattern = "_vs_")) %>%
#   mutate(pair2 = map(pair2, function(x){factor(x, levels = Order.material3) %>% sort() %>% as.character()})) %>%
#   mutate(a = map_chr(pair2, ~.[1]), b = map_chr(pair2, ~.[2])) %>%
#   mutate(a = factor(a, levels = Order.material3), b = factor(b, levels = Order.material3)) %>%
#   arrange(a, b) %>%
#   ggplot(aes(x = a, y = b, fill = R2_mean)) +
#   geom_tile() +
#   scale_fill_gradient(low = "grey90", high = "red", limits = c(0,0.3)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 



DF4 <-
  DF3 %>%
  mutate(data4 = map(data3, 
                     function(x){
                       x %>% 
                         filter(Project == "polymer") %>%
                         separate(Description, into = c("test", "description"), sep = ";") %>%
                         separate(description, into = c("material", "time", "bs")) %>%
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
                         ungroup() %>%
                         dplyr::select(-c(dir, sample, Project, time, time_num)) %>%
                         filter(material2 != "water") %>%
                         dplyr::select(test, material, material2, day, bs, everything()) %>%
                         filter(material != "water") %>%
                         filter(!(test == "test4" | test == "test5" | test == "test5S")) %>%
                         mutate(test = ifelse(test == "test6", "Test.A",
                                              ifelse(test == "test8", "Test.B",
                                                     ifelse(test == "test9", "Test.C",
                                                            ifelse(test == "test10", "Test.D",
                                                                   ifelse(test == "test12", "Test.E", "Other"
                                                                   ))))))})) %>%
  mutate(data5 = map(data4, 
                     function(x){
                       x %>%
                         gather(-test, -material, -material2, -material3, -day, -bs, key = "taxa", value = "prop")
                     })) %>%
  mutate(color.list = map(data5,
                          function(x){
                            tmp.1 <-
                              x %>%
                              nest(-taxa) %>% 
                              mutate(max.prop = map_dbl(data, function(x){x %>% arrange(desc(prop)) %>% dplyr::select(prop) %>% .[[1,1]]}))  %>%
                              arrange(max.prop) %>%
                              dplyr::select(taxa, max.prop)
                            tmp.1.1 <-
                              tmp.1 %>%
                              filter(taxa != "Unassigned")
                            tmp.1.2 <-
                              tmp.1 %>%
                              filter(taxa == "Unassigned")
                            tmp.1.3 <-
                              bind_rows(tmp.1.2, tmp.1.1) 
                            
                            palette <-
                              hue_pal(c = 50)(nrow(tmp.1)) %>%
                              as_tibble() %>%
                              sample_n(nrow(tmp.1)) %>%
                              .[[1]]
                            palette_top12 <-
                              hue_pal()(12)
                            
                            palette[1] <- "grey"
                            for(col in 1:length(palette_top12)){
                              palette[nrow(tmp.1) - col + 1] <- palette_top12[col] 
                            }
                            
                            
                            res <-
                              tmp.1.3 %>%
                              bind_cols(., tibble(color = palette))
                            
                            return(res)
                          })) %>%
  mutate(res.nmds = map(data4,
                        function(x){
                          x %>%
                            dplyr::select(-c(test, material, material2, material3, day, bs)) %>%
                            metaMDS(., k = 2, trymax = 20)
                        })) %>%
  mutate(permanova = map(data4,
                         function(x){
                           tic()
                           tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
                           tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day, bs))
                           res <- adonis2(tmp1 ~ test + material3 + day + bs, data = tmp2, method = "bray")
                           toc()
                           return(res)
                         })) %>%
 mutate(permanova.2 = map(data4,
                         function(x){
                           tic()
                           tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs))
                           tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day, bs))
                           res <- adonis2(tmp1 ~ test * material3 * day * bs, data = tmp2, method = "bray")
                           toc()
                           return(res)
                         })) %>%
#  mutate(pairwise.permanova = map(data4,
#                                  function(x){
#                                    tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#                                    tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day, bs))
#                                    res <- pairwise.adonis(tmp1,tmp2$material2)
#                                    return(res)
#                                  })) %>%
#  mutate(pairwise.permanova.2 = map(data4,
#                                  function(x){
#                                    tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#                                    tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, day, bs))
#                                    res <- pairwise.adonis2(tmp1 ~ material2 * test,data = tmp2)
#                                    return(res)
#                                  })) %>%
  mutate(plot.nmds = map2(data4, res.nmds,
                          function(x, y){
                            x %>%
                              bind_cols(., as_tibble(y$points)) %>%
                              mutate(material2 = factor(material2, levels = Order.material2)) %>%
                              ggplot(aes(x = MDS1, y = MDS2, fill = material3, color = material3, size = day, shape = test)) +
                              geom_point(alpha = 0.5) +
                              scale_shape_manual(values = c("Test.A" = 21, "Test.B" = 22, "Test.C" = 23, "Test.D" = 24, "Test.E" = 25)) +
                              facet_wrap(test ~ .) +
                              theme_bw() +
                              theme(legend.direction = "vertical", legend.box = "vertical") +
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
  mutate(plot.nmds.density = map2(data4, res.nmds,
                          function(x, y){
                            x %>%
                              bind_cols(., as_tibble(y$points)) %>%
                              mutate(material2 = factor(material3, levels = Order.material3)) %>%
                              ggplot(aes(x = MDS1, y = MDS2, color = material3)) + 
                              stat_density_2d(contour_var = "ndensity", binwidth = 0.2) + 
                              facet_wrap(. ~ test) + 
                              theme_bw() +
                              facet_grid(test ~ .) +
                              theme_bw() +
                              theme(legend.direction = "vertical", legend.box = "horizontal") +
                              labs(title = str_c("Stress = ", round(y$stress, 3)),
                                   color = "Material") +
                              scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                            "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                            "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                            "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                            "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                            "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))
                          })) %>%
  mutate(plot_sp.score = map2(res.nmds, color.list,
                              function(x, y){
                                tmp <-
                                  x$species %>% 
                                  as_tibble(rownames = "taxa") %>% 
                                  left_join(., y, by = "taxa") 
                                plot <-
                                  tmp %>%
                                  ggplot(aes(x = MDS1, y = MDS2, color = taxa, size = max.prop)) + 
                                  geom_point() + 
                                  theme_bw() +
                                  theme(legend.position = "none") +
                                  scale_color_manual(values = tmp$color)
                              })) %>%
  mutate(plot_ts = map2(data4, color.list,
                        function(x, y){
                          tmp1 <-
                            x %>%
                            mutate(id = row_number()) %>%
                            mutate(time_unit = ifelse(day < 3, "h", "d")) %>%
                            mutate(time.tmp1 = ifelse(time_unit == "h", round(24 * day, 0), day)) %>%
                            mutate(time.tmp2 = ifelse(time.tmp1 < 10, str_c("0", time.tmp1), time.tmp1)) %>%
                            mutate(time_fct = str_c(time.tmp2, time_unit))
                          time_order <-
                            tmp1 %>%
                            arrange(desc(time_unit), time.tmp2) %>% 
                            .$time_fct %>% 
                            unique()
                          tmp2 <-
                            tmp1 %>%
                            mutate(time_fct = factor(time_fct, levels = time_order)) %>%
                            dplyr::select(-time_unit, -time.tmp1, -time.tmp2) %>%
                            gather(-c(test, material, material2, material3, day, bs, id, time_fct), key = "taxa", value = "prop") %>%
                            mutate(taxa = factor(taxa, levels = y$taxa)) %>%
                            mutate(material3 = factor(material3, levels = Order.material3)) %>%
                            mutate(bs2 = ifelse(bs == "bio", "Biofilm", ifelse(bs == "surf", "Surface", "Not.separated"))) %>%
                            ggplot(aes(x = time_fct, y = prop * 100, fill = taxa)) +
                            geom_bar(stat = "identity", position = "fill") +
                            theme_bw() +
                            theme(legend.position = "none") +
                            scale_fill_manual(values = y$color)+
                            facet_grid(material3 ~ test + bs2 , scales = "free") +
                            #facet_grid2(material2 ~ test, scales = "free", independent = "x") +
                            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                            scale_y_continuous(labels = percent) +
                            ylab("Proportion (%)") +
                            xlab("Incubation Time")
                          
                        })) %>%
  mutate(legend_ts = map(plot_ts,
                         function(x){
                           library(cowplot)
                           tmp <-
                             x$data %>% 
                             nest(-taxa) %>% 
                             mutate(max = map_dbl(data, function(x){max(x$prop)})) %>% 
                             arrange(desc(max)) %>% 
                             filter(taxa != "Unassigned")  %>% 
                             .[1:12,] %>% 
                             mutate(taxa2 = str_split(taxa, pattern = ";")) %>% 
                             mutate(len = map_int(taxa2, ~length(.))) %>% 
                             mutate(taxa3 = map2_chr(taxa2, len, function(x,y){x[y] %>% str_replace(" ", "")})) 
                           if(nrow(tmp) == nrow(tmp %>% nest(-taxa3))){
                             tmp2 <-
                               tmp %>% 
                               mutate(Taxa = factor(taxa3, levels = .$taxa3)) %>% 
                               ggplot(aes(x = Taxa, y = max, fill = Taxa)) + 
                               geom_bar(stat = "identity") +
                               labs(color = "Taxa")
                           } else{
                             tmp2 <- NA
                           }
                           
                           leg.tmp <- get_legend(tmp2)
                           leg <- as_ggplot(leg.tmp)
                           return(leg)
                         }))  %>%
  mutate(plot_ts2 = map2(plot_ts, legend_ts, function(x,y){x + y + plot_layout(ncol = 2, width = c(6, 1))})) %>%
  mutate(res.nmds.3d = map(data4,
                           function(x){
                             set.seed(1234)
                             x %>%
                               dplyr::select(-c(test, material, material2, material3, day, bs)) %>%
                               metaMDS(., k = 3, trymax = 20)
                           })) %>%
  mutate(plot.nmds.3d = map2(data4, res.nmds.3d,
                             function(x, y){
                               tmp <-
                                 x %>%
                                 bind_cols(., as_tibble(y$points)) %>%
                                 mutate(material2 = factor(material2, levels = Order.material2)) %>%
                                 mutate(material3 = factor(material3, levels = Order.material3))
                               pal <- c(rainbow_hcl(n=20 , c= 150, l = 40)[1],
                                        rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                        rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                        rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                        rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                        rainbow_hcl(n=20 , c= 50, l = 70)[16])
                               shape <- c("circle", "square", "diamond", "x", "cross")
                               plot <-
                                 plot_ly(data = tmp,
                                         x=~MDS1, y=~MDS2, z=~MDS3, 
                                         type="scatter3d", mode="markers", 
                                         color=~material3, size = ~day, colors = pal, symbol = ~test, symbols = shape) %>%
                                 layout(title = str_c("Stress = ", round(y$stress,3)))
                             }))  %>%
  mutate(plot.nmds.3d.ts = map2(data4, res.nmds.3d,
                             function(x, y){
                               tmp <-
                                 x %>%
                                 bind_cols(., as_tibble(y$points)) %>%
                                 mutate(material2 = factor(material2, levels = Order.material2)) %>%
                                 mutate(material3 = factor(material3, levels = Order.material3))
                               pal <- c(rainbow_hcl(n=20 , c= 150, l = 40)[1],
                                        rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                        rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                        rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                        rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                        rainbow_hcl(n=20 , c= 50, l = 70)[16])
                               tmp %>%
                                 gather(MDS1, MDS2, MDS3, key = "axis", value = "val") %>%
                                 ggplot(aes(x = day, y = val, color = material3, fill = material3, shape = test)) +
                                 geom_point(alpha = 0.5) +
                                 geom_smooth(method = "lm", se = F) +
                                 theme_bw() +
                                 facet_grid(axis ~ .) +
                                 scale_shape_manual(values = c("Test.A" = 21, "Test.B" = 22, "Test.C" = 23, "Test.D" = 24, "Test.E" = 25)) +
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
                                                              "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) + 
                                 labs(color = "Material", shape = "Test", fill = "Material")
                                 
                             })) %>%
  mutate(plot_ts_each = map(data4,
                            function(x){
                              x %>%
                                gather(-c(test, material, material2, material3, day, bs), key = "taxa", value = prop) %>%
                                nest(-taxa) %>%
                                mutate(median.prop = map_dbl(data, function(x){median(x$prop)})) %>%
                                arrange(desc(median.prop)) %>%
                                mutate(plot = map2(data, taxa, 
                                                   function(xx, yy){
                                                     if(yy == "Unassigned"){
                                                       lab_tmp <- "Unassigned"
                                                     } else {
                                                       lab_tmp <- 
                                                         str_split(yy, pattern = "; ") %>% .[[1]] %>% .[length(.)]
                                                     }
                                                     
                                                     xx %>%
                                                       mutate(material2 = factor(material2, levels = Order.material2)) %>%
                                                       mutate(id = str_c(test, material, bs)) %>%
                                                       ggplot(aes(x = day, y = prop * 100, color = material3, group = id)) +
                                                       geom_point() +
                                                       geom_smooth(method = "lm", se =F) +
                                                       theme_bw() +
                                                       labs(x = "Incubation Day", y = "Proportion (%)", color = "Material") +
                                                       ggtitle(lab_tmp) +
                                                       scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                                     "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                                     "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                                     "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                                     "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                                     "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))
                                                   })) 
                            })) %>%
  mutate(data.bray.distance = map(data4,
                                  function(x){
                                    tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, day, bs))
                                    tmp2 <- x %>% dplyr::select(c(test, material3, day, bs)) %>% mutate(id = row_number() %>% as.character())
                                    tmp3 <- tmp1 %>% vegdist(method = "bray")
                                    tmp4 <-
                                      tmp3 %>%
                                      as.matrix() %>% as_tibble() %>%
                                      bind_cols(tmp2, .) %>%
                                      gather(-c(test,  material3, day, bs, id), key = "id2", value = "distance") %>%
                                      left_join(., tmp2 %>% set_names(c("test.to", "material.to", "day.to", "bs.to", "id2"))) %>%
                                      filter(test == test.to) %>%
                                      filter(id != id2) #%>%
                                    #mutate(pair = map2(material3, material.to, function(x,y){c(as.character(x),as.character(y))%>% sort()})) %>%
                                    #mutate(pair.1 = map_chr(pair, ~.[1]), pair.2 = map_chr(pair, ~.[2])) %>%
                                    #mutate(pair = str_c(pair.1, pair.2, sep = "_vs_"))
                                    
                                  })) %>%
  mutate(plot.bray.distance = map(data.bray.distance,
                                  function(x){
                                    tmp1 <-
                                      x %>%
                                      nest(-material.to, -test) %>%
                                      mutate(multcomp = map(data,
                                                            function(x){
                                                              #tic()
                                                              res <- nparcomp(distance ~ material3, data=x, asy.method = "mult.t",
                                                                              type = "Tukey",alternative = "two.sided",info = FALSE)
                                                              #toc()
                                                              return(res)
                                                            })) %>%
                                      mutate(letters = map(multcomp, func_letter_nparcomp))
                                    
                                    data.letters <- 
                                      tmp1 %>%
                                      dplyr::select(test, material.to, letters) %>%
                                      unnest() %>%
                                      rename(material3 = "param") %>%
                                      mutate(material3 = factor(material3, levels = Order.material3))
                                    
                                    data.distance <-
                                      tmp1 %>%
                                      dplyr::select(test, material.to, data) %>%
                                      unnest() %>%
                                      mutate(material3 = factor(material3, levels = Order.material3))
                                    
                                    ggplot() +
                                      geom_jitter(data = data.distance, aes(x = material3, y = distance, color = material3), alpha = 0.05) +
                                      geom_boxplot(data = data.distance, aes(x = material3, y = distance), fill = NA, outlier.colour = NA) +
                                      geom_text(data = data.letters, aes(x = material3, y = 1.15, label = fin.letter)) +
                                      theme_bw() +  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
                                      facet_grid(test ~ material.to) +
                                      ylim(0, 1.2) +
                                      ylab("Bray-Curtis Index") + xlab("Material") +
                                      scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                    "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                    "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                    "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                    "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                    "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                      theme(legend.position = "none")
                                    
                                    
                                  }))


write_rds(DF4 %>% dplyr::select(category, data4, res.nmds, res.nmds.3d, permanova, permanova.2),
          str_c(Dir_res_microbiome, "/df_raw.rds"))


bray.distance.comparison.raw.rarefied <-
  bind_rows(DF4 %>% filter(category == "asv") %>%
            mutate(type = "raw") %>%
            dplyr::select(type, data.bray.distance),
          DF3_rarefied[1,] %>%
            mutate(type = "rarefied") %>%
            dplyr::select(type, data.bray.distance)) %>%
  unnest() %>%
  mutate(type = factor(type, levels = c("raw", "rarefied"))) %>%
  nest(-type, -material3) %>%
  arrange(material3, type) %>%
  mutate(type2 = str_c(material3, type, sep = "_"))  %>%
  mutate(type2 = factor(type2,levels = .$type2)) %>% unnest() %>%
  nest() %>%
  mutate(plot.bray.distance = map(data,
                                  function(x){
                                    type2.order <-  x$type2 %>% unique() %>% sort()
                                    
                                    tmp1 <-
                                      x %>%
                                      nest(-material.to, -test) %>%
                                      mutate(multcomp = map(data,
                                                            function(x){
                                                              #tic()
                                                              res <- my_own_nparcomp(distance ~ type2, data=x, asy.method = "mult.t",
                                                                              type = "Tukey",alternative = "two.sided",info = FALSE)
                                                              #toc()
                                                              return(res)
                                                            })) %>%
                                      mutate(letters = map(multcomp, func_letter_nparcomp))
                                    
                                    
                                    data.letters <- 
                                      tmp1 %>%
                                      dplyr::select(test, material.to, letters) %>%
                                      unnest() %>%
                                      rename(type2 = "param") %>%
                                      mutate(type2 = factor(type2, levels = type2.order))
                                    
                                    data.distance <-
                                      tmp1 %>%
                                      dplyr::select(test, material.to, data) %>%
                                      unnest() %>%
                                      mutate(material3 = factor(material3, levels = Order.material3),
                                             type2 = factor(type2, levels = type2.order))
                                    
                                    ggplot() +
                                      geom_jitter(data = data.distance, aes(x = type2, y = distance, color = material3), alpha = 0.05) +
                                      geom_boxplot(data = data.distance, aes(x = type2, y = distance), fill = NA, outlier.colour = NA) +
                                      geom_text(data = data.letters, aes(x = type2, y = 1.15, label = fin.letter), size = 2.5) +
                                      theme_bw() +  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
                                      facet_grid(test ~ material.to) +
                                      ylim(0, 1.2) +
                                      ylab("Bray-Curtis Index") + xlab("Material") +
                                      scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                    "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                    "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                    "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                    "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                    "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                      theme(legend.position = "none")
                                    
                                    
                                  }))
ggsave(plot = bray.distance.comparison.raw.rarefied %>% dplyr::select(plot.bray.distance) %>% .[[1,1]] %>% .[[1]],
       file = str_c(Dir_res_microbiome, "/bray_raw.rarefied.png"), height = 6, width = 12) 


# res <- list()
# for(i in 1:nrow(DF4)){
#   tic()
#   tmp <- DF4[i,] %>% dplyr::select(data4) %>% unnest()
#   tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#   tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
#   tmp3 <- pairwise.adonis2(tmp1 ~ material3 + test + day + bs, data = tmp2)
#   res <- append(res, list(tmp3))
#   toc()
# }
# res2 <- list()
# for(i in 1:nrow(DF4)){
#   tic()
#   tmp <- DF4[i,] %>% dplyr::select(data4) %>% unnest()
#   tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#   tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
#   tmp3 <- pairwise.adonis2(tmp1 ~ material3 * test * day * bs, data = tmp2)
#   res2 <- append(res2, list(tmp3))
#   toc()
# }
# 
# 
# DF4_2 <-
#   bind_cols(DF4, tibble(pairwise.permanova = res))  %>%
#   bind_cols(., tibble(pairwise.permanova.2 = res2)) %>%
#   mutate(permanova = map(permanova, function(x)(as_tibble(x, rownames = "factor"))))  %>%
#   mutate(permanova.2 = map(permanova.2, function(x)(as_tibble(x, rownames = "factor"))))  %>%
#   mutate(pairwise.permanova = map(pairwise.permanova, function(x){tibble(pair = names(x), res = x)})) %>%
#   mutate(pairwise.permanova.2 = map(pairwise.permanova.2, function(x){tibble(pair = names(x), res = x)})) 
# 
# save.tmp.2 <-
#   DF4_2 %>%
#   dplyr::select(category, permanova, permanova.2, pairwise.permanova, pairwise.permanova.2)
# 
# write_rds(save.tmp.2, "raw_permanova.rds")
# 
# 
# save.tmp %>%
#   mutate(category = "rarefaction") %>%
#   bind_rows(save.tmp.2 %>% filter(category == "asv")) %>%
#   dplyr::select(category, permanova.2) %>%
#   unnest() %>%
#   filter(factor != "Total") %>%
#   nest(-factor) %>%
#   mutate(factor = str_replace(factor, "3", "")) %>%
#   mutate(factor = str_replace(factor, ":", " x ")) %>%
#   mutate(factor = str_replace(factor, ":", " x ")) %>%
#   mutate(factor = str_replace(factor, ":", " x ")) %>%
#   mutate(factor= factor(factor, levels = .$factor)) %>%
#   unnest() %>%
#   ggplot(aes(x = factor, y = R2, color = category)) +
#   geom_boxplot() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# 
# 
# 
# 
# res3 <- tibble()
# for(i in 1:nrow(DF3_rarefied)){
#   TMP <- DF3_rarefied[i,] %>% dplyr::select(data.rarefied.2) %>% unnest() %>% nest(-test) %>% arrange(test)
#   for(j in 1:nrow(TMP)){
#     tic()
#     tmp <- TMP[j,] %>% unnest()
#     tmp.test <- TMP[j,] %>% dplyr::select(test) %>% .[[1]]
#     tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs))
#     tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs))
#     tmp3 <- try(pairwise.adonis2(tmp1 ~ material3 * day * bs , data = tmp2))
#     tmp4 <- try(pairwise.adonis2(tmp1 ~ material3 * day, data = tmp2))
#     tmp3 <- ifelse(tmp3 == "try-error", "NA", tmp3)
#     tmp4 <- ifelse(tmp4 == "try-error", "NA", tmp4)
#     res.tmp <- tibble(rarefaction = i, test = tmp.test, pairwise.permanova.1 = list(tmp3), pairwise.permanova.2 = list(tmp4))
#     res3 <- bind_rows(res3, res.tmp)
#     rm(tmp, tmp.test, tmp1, tmp2, tmp3, tmp4, res.tmp)
#     toc()
#   }
# }
# 
# 
# 
# 
# perm_asv_each.test <-
#   DF4[[8,3]][[1]] %>%
#   nest(-test) 
# res4 <- list()
# for(i in 1:nrow(perm_asv_each.test)){
#   tic()
#   tmp <- perm_asv_each.test[i,] %>% unnest()
#   tmp.test <- perm_asv_each.test[i,] %>% dplyr::select(test) %>% .[[1]]
#   tmp1 <- tmp %>% dplyr::select(-c(test, material, material2, material3, day, bs)) 
#   tmp2 <- tmp %>% dplyr::select(c(test, material, material2, material3, day, bs)) 
#   tmp3 <- try(pairwise.adonis2(tmp1 ~ material3 * day * bs , data = tmp2))
#   tmp4 <- try(pairwise.adonis2(tmp1 ~ material3 * day, data = tmp2))
#   tmp3 <- ifelse(tmp3 == "try-error", "NA", tmp3)
#   tmp4 <- ifelse(tmp4 == "try-error", "NA", tmp4)
#   res.tmp <- tibble(test = tmp.test, pairwise.permanova.1 = list(tmp3), pairwise.permanova.2 = list(tmp4))
#   res4 <- bind_rows(res4, res.tmp)
#   toc()
# }
# 
# 
# bind_rows(res3, res4)


for(i in 1:nrow(DF4)){
  category <- str_c(i, DF4[i,] %>% dplyr::select(category) %>% .[[1]], sep = ".")
  plot <- DF4[i,] %>% dplyr::select(plot.nmds) %>% .[[1]] %>% .[[1]]
  plot_ts <- DF4[i,] %>% dplyr::select(plot_ts) %>% .[[1]] %>% .[[1]]
  plot_bray <- DF4[i,] %>% dplyr::select(plot.bray.distance) %>% .[[1]] %>% .[[1]]
  ggsave(plot = plot,
         file = str_c(Dir_res_microbiome, "/", category, "_nMDS.png"),
         height = 6, width = 10)
  ggsave(plot = plot + facet_grid() + theme(legend.direction = "vertical", legend.box = "horizontal") ,
         file = str_c(Dir_res_microbiome, "/", category, "_nMDS.2.png"),
         height = 4, width = 7)
  ggsave(plot = plot_ts,
         file = str_c(Dir_res_microbiome, "/", category, "_ts.png"),
         height = 8, width = 16)
  ggsave(plot = plot_bray,
         file = str_c(Dir_res_microbiome, "/", category, "_bray.distance.png"),
         height = 6, width = 16)
}

# for(i in 1:nrow(DF4)){
#   category <- str_c(DF4[i,] %>% mutate(num = row_number()) %>% dplyr::select(num) %>% .[[1]], DF4[i,] %>% dplyr::select(category) %>% .[[1]], sep = ".")
#   tmp <- DF4[i,] %>% dplyr::select(plot_ts_each) %>% .[[1]] %>% .[[1]]
#   dir_tmp <- str_c(Dir_res_microbiome, category, sep = "/")
#   if(!file.exists(dir_tmp)){dir.create(dir_tmp)}
#   
#   for(ii in 1:nrow(tmp)){
#     plot_ts <- tmp[ii,] %>% dplyr::select(plot) %>% .[[1]] %>% .[[1]]
#     taxa_tmp <- tmp[ii,] %>% dplyr::select(taxa) %>% .[[1]]
#     ggsave(plot = plot_ts,
#            file = str_c(dir_tmp, "/", ii, "_",  taxa_tmp, ".png"),
#            height = 3, width = 4)
#   }}


ASV_GLMM <-
  DF4 %>%
  filter(category == "asv") %>%
  dplyr::select(plot.nmds.3d.ts) %>%
  .[[1]] %>% .[[1]] %>%
  .$data %>%
  mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
  mutate(id = str_c(test, material, bs)) %>%
  dplyr::select(test, material3, bs, day, id, axis,val) %>%
  nest(-axis) %>%
  mutate(glmm_mds = map(data,
                        function(x){
                          brm(formula = val ~ day * material3 + (1|id),
                              family = gaussian(link = "identity"),
                              data = x,
                              seed = 1,
                              iter = 3000,
                              warmup = 1000)
                        })) %>%
  mutate(plot = map(glmm_mds,
                    function(x){
                      marginal_effects(x)
                    })) %>%
  mutate(plot_glmm = map2(plot, axis,
                         function(x, y){
                           plot_day <-
                             x$`day` %>%
                             ggplot(aes(x = day, y = estimate__)) +
                             geom_line( size = 1) +
                             geom_ribbon(aes(ymin = lower__, ymax = upper__), color = NA, alpha = 0.1) +
                             theme_bw() +
                             labs(x = "Day", y = y, title = "Day")
                           
                           plot_material <-
                             x$`material3` %>%
                             ggplot(aes(x = material3, y = estimate__, color = material3)) +
                             geom_point() +
                             geom_errorbar(aes(ymax = upper__, ymin = lower__), width = .4) +
                             theme_bw() +
                             labs(x = "Material", y=y, color = "Material", fill = "Material", title = "Material") +
                             scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                           "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                           "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                           "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                           "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                           "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                             theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
                           
                           plot_interaction <-
                             x$`day:material3` %>%
                             ggplot(aes(x = day, y = estimate__, color = material3)) +
                             geom_smooth( se = F) +
                             geom_ribbon(aes(ymin = lower__, ymax = upper__,  fill = material3), color = NA, alpha = 0.1) +
                             theme_bw() +
                             labs(x = "Day", y=y, color = "Material", fill = "Material", title = "Interaction") +
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
                           
                           
                           res <- list(plot_day, plot_material, plot_interaction)
                           #res <- plot_day/plot_material/plot_interaction
                           #res <- plot_interaction
                           return(res)
                           
                         })) %>%
  mutate(plot.glmm.2 = pmap(list(data, plot_glmm, axis),
                            function(x,y,z){
                              ggplot() +
                                geom_point(data = x, aes(x = day, y = val, color = material3, size = day, fill = material3, shape = test), alpha = 0.3) +
                                scale_shape_manual(values = c("Test.A" = 21, "Test.B" = 22, "Test.C" = 23, "Test.D" = 24, "Test.E" = 25)) +
                                geom_line(data = y[[3]]$data, aes(x = day, y = estimate__, color = material3)) +
                                geom_ribbon(data = y[[3]]$data, aes(x = day, y = estimate__, color = material3, ymin = lower__, ymax = upper__,  fill = material3), alpha = 0.1, color = NA) +
                                theme_bw() +
                                labs(x = "Day", y = z, color = "Material", fill = "Material", title = z) +
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
                              
                            }))


plot_asv_glmm <-
  (ASV_GLMM[[1,5]][[1]][[1]] + ASV_GLMM[[1,5]][[1]][[2]] + ASV_GLMM[[1,5]][[1]][[3]])/
  (ASV_GLMM[[2,5]][[1]][[1]] + ASV_GLMM[[2,5]][[1]][[2]] + ASV_GLMM[[2,5]][[1]][[3]])/
  (ASV_GLMM[[3,5]][[1]][[1]] + ASV_GLMM[[3,5]][[1]][[2]] + ASV_GLMM[[3,5]][[1]][[3]])

ggsave(plot = plot_asv_glmm, file = str_c(Dir_res_microbiome, "/result_glmm_asv.pdf"),
       height = 8, width = 12)

Microbiome_output <-
  DF4 %>% 
  dplyr::select(category, data4, res.nmds, permanova, permanova.2)
write_rds(Microbiome_output, str_c(Dir_res_integral, "/Microbiome_output.rds", sep = ""))
# write_rds(DF4, "DF4.rds")

Boxplot <-
  read_rds(str_c(Dir_res_integral, "/Microbiome_output.rds", sep = "")) %>%
  filter(category == "g") %>%
  dplyr::select(data4) %>%
  unnest() %>% 
  mutate(is.PHBH = ifelse(material2 == "PHBH", "PHBH", "Other")) %>% 
  mutate(material3 = factor(material3, levels = Order.material3)) %>%
  gather(-c(test, material, material2, material3, day, bs, is.PHBH), key = "taxa", value = "prop") %>% 
  nest(-taxa) %>%
  mutate(max.prop = map_dbl(data, function(x){x %>% summarize(max(prop)) %>% .[[1]]})) %>%
  filter(max.prop > 0.01) %>%
  mutate(family = str_split(taxa, pattern = "f__"),
         family = map_chr(family, function(x){x[2]}),
         family = str_split(family, pattern = "; g__"),
         family = map_chr(family, function(x){x[1]}),
         genus = str_split(taxa, pattern = "g__"),
         genus = map_chr(genus, function(x){x[2]})) %>%
  mutate(t.test = map(data,
                      function(x){
                        t.test(x %>% filter(is.PHBH == "PHBH") %>% .$prop,
                               x %>% filter(is.PHBH == "Other") %>% .$prop,
                               var.equal = T)})) %>%
  mutate(p.val = map_dbl(t.test, function(x){x$p.value})) %>%
  mutate(anova = map(data,
                     function(x){
                       summary(aov(prop ~ material3, data = x))
                     })) %>%
  mutate(anova.p.val = map_dbl(anova,
                           function(x){
                             x[[1]] %>% .[1,5]
                           })) %>%
  mutate(anova_2way = map(data,
                     function(x){
                       summary(aov(prop ~ test * material3, data = x))
                     })) %>%
  mutate(anova_2way_p.val_test = map_dbl(anova_2way,
                               function(x){
                                 x[[1]] %>% .[1,5]
                               })) %>%
  mutate(anova_2way_p.val_material = map_dbl(anova_2way,
                                         function(x){
                                           x[[1]] %>% .[2,5]
                                         })) %>%
  mutate(anova_2way_p.val_interaction = map_dbl(anova_2way,
                                         function(x){
                                           x[[1]] %>% .[3,5]
                                         })) %>%
  mutate(anova_4way = map(data,
                          function(x){
                            summary(aov(prop ~ test * material3 * day * bs, data = x))
                          })) %>%
  mutate(anova_4way_p_material = map_dbl(anova_4way,
                                         function(x){
                                           x %>%
                                             .[[1]] %>% 
                                             as_tibble(rownames = "factor") %>% 
                                             .[2,6] %>% .[[1]]
                                         })) %>%
  mutate(anova_4way_res = map_chr(anova_4way,
                          function(x){
                            tmp1 <- x[[1]] %>% rownames() %>% as_tibble()
                            tmp2 <- x[[1]] %>%  as_tibble()
                            tmp3 <-
                              bind_cols(tmp1, tmp2) %>%
                              set_names(c("factor", "df", "Sum_seq", "mean_seq", "F.value", "p.value")) %>%
                              mutate(factor = str_replace(factor, pattern = "test", replacement = "Test"),
                                     factor = str_replace(factor, pattern = "material3", replacement = "Material"),
                                     factor = str_replace(factor, pattern = "day", replacement = "Day"),
                                     factor = str_replace(factor, pattern = "bs", replacement = "Biofilm"),
                                     factor = str_trim(factor))
                            tmp4 <- 
                              tmp3 %>%
                              mutate(factor.num = str_split(factor, pattern = ":"),
                                     factor.num = map_int(factor.num, ~length(.))) %>%
                              nest(-factor.num) %>%
                              mutate(tag = map_chr(data,
                                                   function(x){
                                                     x %>%
                                                       filter(p.value < 0.05) %>%
                                                       mutate(factor = str_c(factor, "*")) %>%
                                                       .$factor %>% 
                                                       str_c(collapse = ",") %>%
                                                       str_c("    ", ., sep = "")
                                                   })) 
                            
                            res <- tmp4$tag %>% str_c(collapse = ("\n"))
                            
                            return(res)
                          })) %>%
  mutate(plot = pmap(list(data, genus, anova.p.val), 
                    function(x, y, z){
                      x %>% 
                        ggplot(aes(x = material2, y = prop * 100, color = material3)) + 
                        geom_boxplot() + 
                        theme_bw() + 
                        ylab("Proportion (%)") + 
                        xlab("Material") +
                        labs(color = "Material") +
                        scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                      "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                      "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                      "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                      "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                      "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                        ggtitle(str_c(y, " ; ANOVA ", ifelse(z < 0.05, "p < 0.05", str_c("p = ", round(z, 2)))))})) %>%
  mutate(plot2 = pmap(list(plot, family, genus, anova_2way_p.val_test, anova_2way_p.val_material, anova_2way_p.val_interaction), 
                     function(XX, YY, ZZ, x, y, z){
                       XX +
                         facet_grid(test ~ .) +
                         ggtitle(str_c(YY, " ", ZZ, "\n Test ", 
                                       ifelse(x < 0.05, "p < 0.05", str_c("p = ", round(x, 2))),
                                       "\n Material ",
                                       ifelse(y < 0.05, "p < 0.05", str_c("p = ", round(y, 2))),
                                       "\n Interaction ",
                                       ifelse(z < 0.05, "p < 0.05", str_c("p = ", round(z, 2)))
                                       ))
                         })) %>%
  mutate(plot3 = pmap(list(plot, family, genus, anova_2way_p.val_test, anova_2way_p.val_material, anova_2way_p.val_interaction), 
                      function(XX, YY, ZZ, x, y, z){
                        XX$data %>%
                          mutate(id = str_c(test, material, bs, sep = "_")) %>%
                          mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
                          ggplot(aes(x = day, y = 100* prop, color = material3, group = id, shape = bs2)) +
                          geom_line(alpha = 0.4) +
                          geom_point(alpha = 0.4) +
                          #geom_smooth(se = F, alpha = 0.4) +
                          facet_grid(test ~ .) +
                          theme_bw() +
                          ylab("Proportion (%)") + 
                          xlab("Incubation Day") +
                          labs(color = "Material", shape = "Biofilm/Surface") +
                          scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                        "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                        "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                        "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                        "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                        "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                          ggtitle(str_c(YY, " ", ZZ))
                      })) %>%
  
  #mutate(plot4 = map2(plot2, plot3, function(x, y){x + theme(legend.position = "none") + y + ggtitle("")})) %>%
  mutate(plot4 = pmap(list(plot2, plot3, family, genus, 
                           anova_2way_p.val_test , anova_2way_p.val_material, anova_2way_p.val_interaction, anova_4way_res), 
                      function(XX, YY, x, y, a, b, c, d){
                        XX + theme(legend.position = "none") + theme(plot.title = element_blank()) + 
                          YY + theme(plot.title = element_blank()) + 
                          plot_annotation(title = str_c(x, " ", y),
                                          subtitle = str_c("4-way ANOVA;\n",d))
                      })) %>%
  arrange(desc(max.prop)) %>%
  mutate(median.PHBH = map_dbl(data, function(x){x %>% filter(is.PHBH == "PHBH") %>% .$prop %>% median()})) %>%
  mutate(median.other = map_dbl(data, function(x){x %>% filter(is.PHBH != "PHBH") %>% .$prop %>% median()})) %>%
  mutate(diff = median.PHBH - median.other) %>%
  arrange(anova.p.val) %>%
  mutate(multicomp = map(data,
                         function(x){
                           tmp1 <- aov(prop ~ material3, data = x)
                           tmp2 <- glht(tmp1,linfct = mcp(material3 = "Tukey"))
                           summary(tmp2)
                         })) %>%
  mutate(res.mult = map(multicomp,
                        function(x){
                          x$test %>% 
                            .[c(3:6)]  %>% 
                            as_tibble() %>% 
                            bind_cols(x$linfct %>% rownames() %>% as_tibble() %>% set_names("comparison"),
                                      .)
                        })) %>%
  mutate(is.only.hl = map(res.mult,
                              function(x){
                                tmp <-
                                  x %>%
                                  filter(pvalues < 0.05) %>% 
                                  separate(comparison, into = c("a", "b")) %>% 
                                  gather(a,b, key = "ab", value = "material") %>% 
                                  group_by(material) %>% 
                                  count(material) %>% 
                                  filter(n == 4)
                                return(tmp)
                              })) %>%
  arrange(desc(max.prop)) %>%
  mutate(genus = ifelse(is.na(genus), "Unassigned", genus)) %>%
  mutate(plot.ts_early = pmap(list(data, family, genus),
                              function(x, y, z){
                                color_5 <- hue_pal()(5)
                                x %>%
                                  mutate(id = str_c(test, material, bs, sep = "_")) %>%
                                  mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
                                  filter(test == "Test.D" | test == "Test.E") %>%
                                  ggplot(aes(x = day, y = 100* prop, color = material3, group = id, shape = bs2)) +
                                  geom_line(alpha = 0.4) +
                                  geom_point(alpha = 0.4) +
                                  #geom_smooth(se = F, alpha = 0.4) +
                                  facet_grid(test ~ .) +
                                  theme_bw() +
                                  ylab("Proportion (%)") + 
                                  xlab("Incubation Day") +
                                  labs(color = "Material", shape = "Biofilm/Surface") +
                                  ggtitle(str_c(y, " ", z)) +
                                  scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                                "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                                "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                                "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                                "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                                "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
                              }))
  


Boxplot2 <- Boxplot %>% filter(anova_4way_p_material < 0.05) 

Dir_res_microbiome_boxplot <- str_c(CWD, "/res/microbiome/boxplot_g")

if(!file.exists(Dir_res_microbiome_boxplot)){
  dir.create(Dir_res_microbiome_boxplot)
}


for (i in 1:nrow(Boxplot2)){
  name.g <-  Boxplot2[i,] %>% dplyr::select(genus) %>% .[[1]]
#  plot <- Boxplot2[i,] %>% dplyr::select(plot2) %>% .[[1]] %>% .[[1]]
  plot2 <- Boxplot2[i,] %>% dplyr::select(plot4) %>% .[[1]] %>% .[[1]]
  plot3 <- Boxplot2[i,] %>% dplyr::select(plot.ts_early) %>% .[[1]] %>% .[[1]] + theme(legend.direction = "vertical", legend.box = "vertical")
  
#  ggsave(plot = plot, file = str_c(Dir_res_microbiome_boxplot, "/boxplot_", i ,"_", name.g, ".png"), height = 8, width = 5)
  ggsave(plot = plot2, file = str_c(Dir_res_microbiome_boxplot, "/combineplot_", i ,"_", name.g, ".png"), height = 6.4, width = 6.4)
  ggsave(plot = plot3, file = str_c(Dir_res_microbiome_boxplot, "/ts.early_", i ,"_", name.g, ".png"), height = 3, width = 6.5)
}

# 
# nMDS_each_test_f <-
#   Microbiome_output %>%
#   filter(category == "f") %>%
#   dplyr::select(data4) %>%
#   unnest() %>%
#   nest(-test) %>% 
#   arrange(test) %>% 
#   mutate(res.nmds = map(data,
#                         function(x){
#                           metaMDS(x %>% dplyr::select(-c(material, material2, material3, day, bs)), k = 2, trymax = 20)
#                         })) %>%
#   mutate(plot.nmds = pmap(list(data, res.nmds, test),
#                           function(x, y, z){
#                             color_5 <- hue_pal()(5)
#                             x %>%
#                               bind_cols(., as_tibble(y$points)) %>%
#                               arrange(day) %>%
#                               mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
#                               ggplot(aes(x = MDS1, y = MDS2, color = material3, size = day, shape = bs2)) +
#                               scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                             "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                             "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                             "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                             "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                             "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                               scale_size_continuous(limits = c(0, 56)) +
#                               geom_point(alpha = 0.5) +
#                               #facet_grid(test ~ .) +
#                               theme_bw() +
#                               guides(shape = guide_legend(order = 3),
#                                      color = guide_legend(order = 1),
#                                      size = guide_legend(order = 2)) +
#                               theme(legend.direction = "vertical", legend.box = "horizontal") +
#                               ggtitle(str_c(z, "; stress = ", y$stress %>% round(.,2), sep = "")) +
#                               labs(size = "Day", color = "Material", shape = "Biofilm/Surface")
#                           })) %>%
#     mutate(pairwise.permanova = map(data,
#                                     function(x){
#                                       tmp1 <- x %>% dplyr::select(-c(material, material2, material3, day, bs)) 
#                                       tmp2 <- x %>% dplyr::select(c(material, material2, material3, day, bs))
#                                       res <- pairwise.adonis(tmp1,tmp2$material3)
#                                       return(res)
#                                     }))
# 
# legend <-
#   nMDS_each_test_f[[1,4]][[1]] %>% get_legend() %>% as_ggplot() 
# 
# plot_nmds_each_f <-
#   nMDS_each_test_f[[1,4]][[1]] + theme(legend.position = "none") +
#   nMDS_each_test_f[[2,4]][[1]] + theme(legend.position = "none") +
#   nMDS_each_test_f[[3,4]][[1]] + theme(legend.position = "none") +
#   nMDS_each_test_f[[4,4]][[1]] + theme(legend.position = "none") +
#   nMDS_each_test_f[[5,4]][[1]] + theme(legend.position = "none") +
#   legend
# 
# ggsave(plot = plot_nmds_each_f ,
#        file = str_c(Dir_res_microbiome, "/", "f_nMDS.each.png"),
#        width = 12, height = 8)


# 
# ### bar plot of top 12 (max proportion) and shannon
# bar <-
#   Microbiome_output[-1,] %>%
#   mutate(bar = map(data4,
#                    function(XX){
#                      XX %>%
#                        gather(-c(test, material, material2, material3, day, bs), key = "taxa", value = "prop") %>%
#                        mutate(material2 = factor(material2, levels = Order.material2)) %>%
#                        mutate(material3 = factor(material3, levels = Order.material3)) %>%
#                        nest(-taxa) %>%
#                        mutate(max = map_dbl(data, function(x){max(x$prop)}),
#                               median = map_dbl(data, function(x){median(x$prop)})) %>%
#                        arrange(desc(max)) %>%
#                        filter(taxa != "Unassigned") %>%
#                        .[1:12,] %>%
#                        mutate(taxa2 = str_split(taxa, pattern = ";")) %>%
#                        mutate(len = map_dbl(taxa2, function(x){length(x)})) %>%
#                        mutate(taxa3 = map2_chr(taxa2, len, function(x, y){str_c(x[y-1], x[y], sep = ";")})) %>%
#                        mutate(taxa = factor(taxa, levels = .$taxa)) %>%
#                        mutate(taxa3 = factor(taxa3, levels = .$taxa3)) %>%
#                        unnest(data) %>%
#                        ggplot(aes(x = material2, y = prop * 100, color = material3)) +
#                        geom_boxplot() +
#                        facet_wrap(taxa3 ~ ., labeller = label_wrap_gen(5), scales = "free") +
#                        #theme(legend.position  = "none") +
#                        theme_bw() +
#                        scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                      "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                      "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                      "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                      "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                      "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                        labs(color = "Material", y = "Proprtion (%)") +
#                        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
#                    })) %>%
#   mutate(bar2 = map(data4,
#                    function(XX){
#                      XX %>%
#                        gather(-c(test, material, material2, material3, day, bs), key = "taxa", value = "prop") %>%
#                        mutate(material2 = factor(material2, levels = Order.material2)) %>%
#                        mutate(material3 = factor(material3, levels = Order.material3)) %>%
#                        nest(-taxa) %>%
#                        mutate(max = map_dbl(data, function(x){max(x$prop)}),
#                               median = map_dbl(data, function(x){median(x$prop)})) %>%
#                        arrange(desc(median)) %>%
#                        filter(taxa != "Unassigned") %>%
#                        .[1:12,] %>%
#                        mutate(taxa2 = str_split(taxa, pattern = ";")) %>%
#                        mutate(len = map_dbl(taxa2, function(x){length(x)})) %>%
#                        mutate(taxa3 = map2_chr(taxa2, len, function(x, y){str_c(x[y-1], x[y], sep = ";")})) %>%
#                        mutate(taxa = factor(taxa, levels = .$taxa)) %>%
#                        mutate(taxa3 = factor(taxa3, levels = .$taxa3)) %>%
#                        unnest(data) %>%
#                        ggplot(aes(x = material2, y = prop * 100, color = material3)) +
#                        geom_boxplot() +
#                        facet_wrap(taxa3 ~ ., labeller = label_wrap_gen(5), scales = "free") +
#                        #theme(legend.position  = "none") +
#                        theme_bw() +
#                        scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                      "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                      "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                      "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                      "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                      "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                        labs(color = "Material", y = "Proprtion (%)") +
#                        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
#                    })) %>%
#   mutate(diversity = map(data4, 
#                          function(x){
#                            tmp1 <-
#                              x %>% 
#                              dplyr::select(-c(test, material, material2, material3, day, bs)) %>%
#                              diversity(index = "shannon", base = 2) %>%
#                              as_tibble() %>%
#                              set_names("Shannon")
#                            tmp2 <-
#                              x %>% 
#                              dplyr::select(-c(test, material, material2, material3, day, bs)) %>%
#                              diversity(index = "simpson", base = 2) %>%
#                              as_tibble() %>%
#                              set_names("Simpson")
#                            tmp3 <- 
#                              x %>%
#                              bind_cols(., tmp1, tmp2) %>%
#                              mutate(material2 = factor(material2 , levels = Order.material2)) %>%
#                              mutate(material3 = factor(material3 , levels = Order.material3))
#                            return(tmp3)
#                          })) %>%
#   mutate(plot.diversity = map(diversity,
#                               function(x){
#                                 x %>%
#                                   mutate(id = str_c(test, material, bs)) %>%
#                                   mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
#                                   gather(Shannon, Simpson, key = "Diversity", value = "val") %>%
#                                   ggplot(aes(x = day, y = val, color = material3, group = id, shape = bs2)) +
#                                   geom_point(alpha = 0.4) +
#                                   #geom_smooth(method = "lm", method.args= list(degree = 2), se = F) +
#                                   geom_line(alpha = 0.4) +
#                                   facet_grid2(test ~ Diversity, scales = "free", independent = "y") +
#                                   scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                                 "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                                 "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                                 "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                                 "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                                 "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                                   theme_bw() +
#                                   labs(y = expression(paste(alpha, "-diversity")), x = "Incubation Day", color = "Material", shape = "Biofilm/Surface")
#                               })) %>%
#   mutate(plot.diversity_selected = map(diversity,
#                                        function(x){
#                                          color_5 <- hue_pal()(5)
#                                          x %>%
#                                            filter(test == "Test.D" | test == "Test.E") %>%
#                                            mutate(id = str_c(test, material, bs)) %>%
#                                            mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
#                                            gather(Shannon, Simpson, key = "Diversity", value = "val") %>%
#                                            ggplot(aes(x = day, y = val, color = material3, group = id, shape = bs2)) +
#                                            geom_point(alpha = 0.4) +
#                                            #geom_smooth(method = "lm", method.args= list(degree = 2), se = F) +
#                                            geom_line(alpha = 0.4) +
#                                            facet_grid2(test ~ Diversity, scales = "free", independent = "y") +
#                                            theme_bw() +
#                                            xlim(0,14) +
#                                            labs(y = expression(paste(alpha, "-diversity")), x = "Incubation Day", color = "Material", shape = "Biofilm/Surface") +
#                                            scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                                          "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                                          "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                                          "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                                          "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                                          "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
#                                        })) 
# 
# 
# for(i in 1:nrow(bar)){
#   category <- bar[i,] %>% dplyr::select(category) %>% .[[1]] 
#   plot <- bar[i,] %>% dplyr::select(bar) %>% .[[1]] %>% .[[1]]
#   plot2 <- bar[i,] %>% dplyr::select(plot.diversity) %>% .[[1]] %>% .[[1]]
#   plot3 <- bar[i,] %>% dplyr::select(plot.diversity_selected) %>% .[[1]] %>% .[[1]]
#   ggsave(plot = plot,
#          file = str_c(Dir_res_microbiome, "/top12_boxplot_", i+1, ".", category, ".png"),
#          width = 8, height = 4)
#   ggsave(plot = plot2,
#          file = str_c(Dir_res_microbiome, "/Diversity_", i+1, ".", category, ".png"),
#          width = 6, height = 4)
#   ggsave(plot = plot3,
#          file = str_c(Dir_res_microbiome, "/Diversity.selected_", i+1, ".", category, ".png"),
#          width = 6, height = 3)
# }


# Target
PHBH.degrader <-
  c("Marisediminitalea", "Aestuariibacter", "Alcanivorax", "Alcaligenes", "Pseudoalteromonas", "Alteromonas", "Bacillus",
    "Comamonas", "Enterobacter", "Aliiglaciecola", "Gracilibacillus", "Marinobacter", "Nocardiopsis", "Pseudoalteromonas",
    "Psychrobacillus", "Rheinheimera", "Shewanella", "Streptomyces") # ref: https://www.nature.com/articles/s41428-020-00396-5

DF_PHBH.degrader <-
  Microbiome_output[[6,2]][[1]] %>%
  gather(-c(test, material, material2, material3, day, bs), key = "taxa", value = "prop") %>%
  nest(-taxa) %>%
  filter(str_detect(taxa, paste(PHBH.degrader, collapse = "|"))) %>%
  mutate(max = map_dbl(data, function(x){max(x$prop)})) %>% 
  arrange(desc(max)) %>%
  mutate(plot = map(data, 
                    function(x){
                      x %>% 
                        mutate(material2 = factor(material2, levels = Order.material2)) %>%
                        mutate(material3 = factor(material3, levels = Order.material3)) %>%
                        ggplot(aes(x = material2, y = prop * 100, color = material3)) + 
                        geom_boxplot() +
                        theme_bw() +
                        labs(y = "Proportion (%)", x = "Material", color = "Material") +
                        scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                      "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                      "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                      "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                      "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                      "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
                        })) 


# 
# 
# ### final figure
# figure.3_1 <-
#   DF4[5,] %>% 
#   dplyr::select(plot.nmds)  %>% 
#   .[[1]] %>% .[[1]] + 
#   facet_grid() +                               
#   theme(legend.position = "none") +
#   labs(x = "nMDS1", y = "nMDS2") 
# 
# figure.3_2.1 <-
#   Family_GLMM[1,] %>% dplyr::select(plot.glmm.2) %>% .[[1]] %>% .[[1]] +
#   scale_x_continuous(breaks = seq(0,60,10)) +
#   scale_y_continuous(breaks = seq(-1.5,1.5,0.5)) +
#   theme(legend.position = "none", 
#         plot.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()) +
#   labs(y = "nMDS1")
# figure.3_2.2 <-
#   Family_GLMM[2,] %>% dplyr::select(plot.glmm.2) %>% .[[1]] %>% .[[1]] +
#   scale_x_continuous(breaks = seq(0,60,10)) +
#   scale_y_continuous(breaks = seq(-1.5,1.5,0.5)) +
#   theme(legend.position = "none", 
#         plot.title = element_blank()) +
#   labs(y = "nMDS2")
# library(ggpubr)
# legend_source <- 
#   Family_GLMM[1,] %>% dplyr::select(plot.glmm.2) %>% .[[1]] %>% .[[1]] +
#   theme(legend.direction = "vertical", legend.box = "horizontal") +
#   labs(color = "Material", fill = "Material", shape = "Test", size = "Day") +
#   guides(shape = guide_legend(order = 1),
#          color = guide_legend(order = 2),
#          size = guide_legend(order = 3),
#          fill = F) 
# legend <- legend_source  %>% get_legend()
# plot_legend<- as_ggplot(legend)
# 
# 
# Figure.3 <- {figure.3_1 + {figure.3_2.1/figure.3_2.2}} + plot_legend
# ggsave(plot = Figure.3,
#        file = str_c(Dir_res_microbiome, "/Figure.3.pdf"),
#        width = 8.5, height = 4)



######################
### Figs for thesis
######################

# save manually because plotly cannot export 3d plot in pdf
Figure.2 <- DF4 %>% filter(category == "asv") %>% dplyr::select(plot.nmds.3d) %>% .[[1]] %>% .[[1]]

Figure.S3 <- DF4 %>% filter(category == "asv") %>% dplyr::select(plot.bray.distance) %>% .[[1]] %>% .[[1]]
Figure.S4 <- plot_asv_glmm  
Figure.S5 <- DF3_rarefied[1,] %>% dplyr::select(plot_nmds) %>% .[[1]] %>% .[[1]]
Figure.S6 <- bray.distance.comparison.raw.rarefied %>% dplyr::select(plot.bray.distance) %>% .[[1,1]] %>% .[[1]]
Figure.S7 <- DF4 %>% filter(category == "p") %>% dplyr::select(plot_ts2) %>% .[[1]] %>% .[[1]]
Figure.S8 <- DF4 %>% filter(category == "c") %>% dplyr::select(plot_ts2) %>% .[[1]] %>% .[[1]]
Figure.S9 <- DF4 %>% filter(category == "o") %>% dplyr::select(plot_ts2) %>% .[[1]] %>% .[[1]]
Figure.S10 <- DF4 %>% filter(category == "f") %>% dplyr::select(plot_ts2) %>% .[[1]] %>% .[[1]]

Figure.S12.pre <-
  Boxplot %>%
  filter(str_detect(taxa, "Thalassolituus")|str_detect(taxa, "Simiduia")|str_detect(taxa, "Pseudomonas")|str_detect(taxa, "Pseudarcobacter")) %>%
  dplyr::select(taxa, plot3) 
Figure.S12 <-
  gridExtra::grid.arrange(
    Figure.S12.pre[[3,2]][[1]] + theme(legend.position = "none"), 
    Figure.S12.pre[[4,2]][[1]] + theme(legend.position = "none"), 
    Figure.S12.pre[[1,2]][[1]] + theme(legend.position = "none"),
    Figure.S12.pre[[2,2]][[1]] + theme(legend.position = "none"),
    Figure.S12.pre[[3,2]][[1]] %>% get_legend() %>% as_ggplot(), 
    layout_matrix = rbind(c(1,1,2,2,5), c(3,3,4,4,5)) )


ggsave(plot = Figure.S3, file = str_c(Dir_res_final_fig, "/Figure.S3.pdf"),
       height = 6, width = 10)
ggsave(plot = Figure.S4, file = str_c(Dir_res_final_fig, "/Figure.S4.pdf"),
       height = 8, width = 12)
ggsave(plot = Figure.S5, file = str_c(Dir_res_final_fig, "/Figure.S5.pdf"),
       height = 4, width = 7)
ggsave(plot = Figure.S6, file = str_c(Dir_res_final_fig, "/Figure.S6.pdf"),
       height = 6, width = 12)
ggsave(plot = Figure.S7, file = str_c(Dir_res_final_fig, "/Figure.S7.pdf"),
       height = 8, width = 20)
ggsave(plot = Figure.S8, file = str_c(Dir_res_final_fig, "/Figure.S8.pdf"),
       height = 8, width = 20)
ggsave(plot = Figure.S9, file = str_c(Dir_res_final_fig, "/Figure.S9.pdf"),
       height = 8, width = 20)
ggsave(plot = Figure.S10, file = str_c(Dir_res_final_fig, "/Figure.S10.pdf"),
       height = 8, width = 20)
ggsave(plot = Figure.S12, file = str_c(Dir_res_final_fig, "/Figure.S12.pdf"),
       height = 8, width = 9)


