rm(list = ls())

library(tidyverse)
library(lubridate)
library(scales)
library(vegan)
library(pairwiseAdonis)
library(jsonlite)
library(ggh4x) # for extended facet in ggplot2
library(tictoc)
library(multcomp)
library(nparcomp)
source("helper_functions.r")

library(rstan)
library(brms)
library(tidybayes)
library(ggdist)
library(ggridges)
library(patchwork)
library(colorspace)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

CWD <- getwd()
Dir_data_miseq <- str_c(CWD, "/data") #"D:/data/Miseq/Qiime2"
Dir_res_picrust2 <- str_c(CWD, "/res/picrust2")
Dir_res_integral <- str_c(CWD, "/res/integral")
Dir_res_final_fig <- str_c(CWD, "/res/final.fig")

if(!file.exists(Dir_res_picrust2)){
  dir.create(Dir_res_picrust2)
}
if(!file.exists(Dir_res_integral)){
  dir.create(Dir_res_integral)
}
if(!file.exists(Dir_res_final_fig)){
  dir.create(Dir_res_final_fig)
}

##################################################
### make KO list

KO_list <-
  fromJSON(str_c(Dir_data_miseq, "/ko00001.json")) %>%
  pluck("children") %>% 
  as_tibble() %>%
  rename(Group.1 = "name") %>%
  mutate(children = map(children,
                        function(x){
                          as_tibble(x) %>% 
                            rename(Group.2 = "name") %>%
                            mutate(children = map(children,
                                                  function(xx){
                                                    as_tibble(xx) %>%
                                                      rename(Group.3 = "name")
                                                  }))})) %>%
  unnest() %>%
  unnest() %>%
  unnest() %>%
  separate(name, into = c("tmp", "description"), sep = "; ") %>%
  separate(tmp, into = c("KO", "gene.name"), sep = "  ")
   


##################################################

Order.material2 <-
  c("PHBH", "PCL", "PBSA", "PBS", "PBAT")
Order.material3 <-
  c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT")

DF_PICRUSt <- 
  read_rds(str_c(Dir_data_miseq, "/DF_PICRUSt.rds", sep = ""))

DF_PICRUSt_2 <-
  DF_PICRUSt %>%
  filter(db == "data_ko.2") %>%
  mutate(data3 = map(data2,
                     function(x){
                       x %>%
                         filter(Project == "polymer")  %>%
                         separate(Description, into = c("test", "info"), sep = ";") %>%
                         filter(test == "test6" | test == "test8" | test == "test9" | test == "test10" | test == "test12") %>%
                         separate(info, into = c("material", "time", "bs"))  %>%
                         mutate(time_num = str_sub(time, 1, -2) %>% as.double()) %>%
                         mutate(day = ifelse(str_detect(time, "h"), time_num /24, 
                                             ifelse(str_detect(time, "d"), time_num,
                                                    ifelse(str_detect(time, "w"), time_num*7, NA)))) %>%
                         dplyr::select(-c(dir, time_num, sample, Project)) %>%
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
                         dplyr::select(test,  material, material2, material3, time,  bs, day, everything()) %>%
                         mutate_at(vars(-c(test,  material, material2, material3, time,  bs, day)), ~ifelse(is.na(.), 0, .)) %>%
                         mutate(test = ifelse(test == "test6", "Test.A",
                                              ifelse(test == "test8", "Test.B",
                                                     ifelse(test == "test9", "Test.C",
                                                            ifelse(test == "test10", "Test.D",
                                                                   ifelse(test == "test12", "Test.E", "Other"
                                                                   ))))))
                     })) %>%
  dplyr::select(-data2) %>%
  mutate(nmds = map(data3,
                    function(x){
                      x %>%
                        dplyr::select(-c(test, material, material2, material3, time, bs, day )) %>%
                        metaMDS(., k = 2, trymax = 20)
                    })) %>%
  mutate(df.nmds = map2(data3, nmds,
                          function(x, y){
                            x %>%
                              dplyr::select(c(test, material, material2, material3, time, bs, day )) %>%
                              bind_cols(., y$points %>% as_tibble()) %>%
                              mutate(id2 = str_c(test, material, bs, sep = "_")) 
                          })) %>%
  mutate(plot.nmds = map2(df.nmds, nmds,
                          function(x, y){
                            x %>%
                              ggplot(aes(x = MDS1, y = MDS2, color = material3, fill = material3, size = day, shape = test)) +
                              geom_point(alpha = 0.5) +
                              scale_shape_manual(values = c("Test.A" = 21, "Test.B" = 22, "Test.C" = 23, "Test.D" = 24, "Test.E" = 25)) +
                              theme_bw() +
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
  mutate(permanova2 = map(data3,
                          function(x){
                            tmp1 <- x %>% dplyr::select(-c(test, material, material2, material3, time, day, bs)) 
                            tmp2 <- x %>% dplyr::select(c(test, material, material2, material3, time, day, bs))
                            res <- adonis2(tmp1 ~ test*material3*day*bs, data = tmp2, method = "bray")
                            return(res)
                          })) %>%
  mutate(glmm_mds1 = map(df.nmds,
                    function(x){
                      brm(formula = MDS1 ~ day * material3 + (1|id2),
                          family = gaussian(link = "identity"),
                          data = x,
                          seed = 1,
                          iter = 3000,
                          warmup = 1000)
                    })) %>%
  mutate(glmm_mds2 = map(df.nmds,
                         function(x){
                           brm(formula = MDS2 ~ day * material3 + (1|id2),
                               family = gaussian(link = "identity"),
                               data = x,
                               seed = 1,
                               iter = 3000,
                               warmup = 1000)
                         })) %>%
  mutate(plot1_mds1 = map(glmm_mds1, function(x){plot(x)})) %>%
  mutate(plot1_mds2 = map(glmm_mds2, function(x){plot(x)})) %>%
  mutate(fixed_effect_mds1 = map(glmm_mds1,
                            function(x){
                              x %>%
                                summary() %>%
                                .$fixed %>%
                                as_tibble(rownames = "factor")
                            })) %>%
  mutate(fixed_effect_mds2 = map(glmm_mds2,
                                 function(x){
                                   x %>%
                                     summary() %>%
                                     .$fixed %>%
                                     as_tibble(rownames = "factor")
                                 })) %>%
  mutate(mcmc_mds1 = map(glmm_mds1,
                     function(x){
                       x %>%
                         spread_draws(b_Intercept, b_day, `b_material3PHBH_HH6%`, b_material3PCL, b_material3PBSA, b_material3PBS, b_material3PBAT, 
                                      `b_day:material3PHBH_HH6%`, `b_day:material3PCL`, `b_day:material3PBSA`, `b_day:material3PBS`, `b_day:material3PBAT`) %>%
                         gather(-c(.chain, .iteration, .draw), key = "factor", value = "val") %>%
                         mutate(factor = str_replace(factor, pattern = "material3", replace = "material")) %>%
                         mutate(factor2 = ifelse(str_detect(factor, "b_day:material"), "Interaction", 
                                                 ifelse(str_detect(factor, "b_material"), "Material", 
                                                        ifelse(factor == "b_day", "Day",
                                                               ifelse(factor == "b_Intercept", "Intercept", "other"))))) %>%
                         mutate(factor2 = factor(factor2, levels = c("Intercept", "Day", "Material", "Interaction"))) %>%
                         mutate(factor = str_replace(factor, pattern = "b_", replacement = "")) %>% 
                         mutate(factor = str_replace(factor, pattern = "material", replacement = "")) %>%
                         mutate(factor = str_replace(factor, pattern = "day", replacement = "Day")) %>%
                         ggplot(aes(x = val, y = factor, fill = factor(stat(quantile)))) +
                         stat_density_ridges(geom = "density_ridges_gradient",
                                             calc_ecdf = TRUE,
                                             quantiles = c(0.025, 0.975)) +
                         scale_fill_manual(name = "Probability", 
                                           values = c("#A0A0A0A0", "#00BFC4", "#A0A0A0A0"),
                                           labels = c("0 - 2.5 %", "2.5 - 97.5 %", "97.5 - 100 %")) +
                         facet_grid2(factor2 ~ ., scales = "free", independent = "x") +
                         theme_classic() +
                         geom_vline(xintercept = 0, color = "#F8766D", linetype = 2, size = 1) +
                         theme( strip.background = element_blank()) +
                         xlab("Estimate") +
                         ylab("Factor")
                     })) %>%
  mutate(mcmc_mds2 = map(glmm_mds2,
                         function(x){
                           x %>%
                             spread_draws(b_Intercept, b_day, `b_material3PHBH_HH6%`, b_material3PCL, b_material3PBSA, b_material3PBS, b_material3PBAT, 
                                          `b_day:material3PHBH_HH6%`, `b_day:material3PCL`, `b_day:material3PBSA`, `b_day:material3PBS`, `b_day:material3PBAT`) %>%
                             gather(-c(.chain, .iteration, .draw), key = "factor", value = "val") %>%
                             mutate(factor = str_replace(factor, pattern = "material3", replace = "material")) %>%
                             mutate(factor2 = ifelse(str_detect(factor, "b_day:material"), "Interaction", 
                                                     ifelse(str_detect(factor, "b_material"), "Material", 
                                                            ifelse(factor == "b_day", "Day",
                                                                   ifelse(factor == "b_Intercept", "Intercept", "other"))))) %>%
                             mutate(factor2 = factor(factor2, levels = c("Intercept", "Day", "Material", "Interaction"))) %>%
                             mutate(factor = str_replace(factor, pattern = "b_", replacement = "")) %>% 
                             mutate(factor = str_replace(factor, pattern = "material", replacement = "")) %>%
                             mutate(factor = str_replace(factor, pattern = "day", replacement = "Day")) %>%
                             ggplot(aes(x = val, y = factor, fill = factor(stat(quantile)))) +
                             stat_density_ridges(geom = "density_ridges_gradient",
                                                 calc_ecdf = TRUE,
                                                 quantiles = c(0.025, 0.975)) +
                             scale_fill_manual(name = "Probability", 
                                               values = c("#A0A0A0A0", "#00BFC4", "#A0A0A0A0"),
                                               labels = c("0 - 2.5 %", "2.5 - 97.5 %", "97.5 - 100 %")) +
                             facet_grid2(factor2 ~ ., scales = "free", independent = "x") +
                             theme_classic() +
                             geom_vline(xintercept = 0, color = "#F8766D", linetype = 2, size = 1) +
                             theme( strip.background = element_blank()) +
                             xlab("Estimate") +
                             ylab("Factor")
                         })) %>%
  mutate(plot2_mds1 = map(glmm_mds1,
                     function(x){
                       marginal_effects(x)
                     })) %>%
  mutate(plot2_mds2 = map(glmm_mds2,
                          function(x){
                            marginal_effects(x)
                          })) %>%
  mutate(plot2_mds1.1 = map(plot2_mds1,
                          function(x){
                            plot_day <-
                              x$`day` %>%
                              ggplot(aes(x = day, y = estimate__)) +
                              geom_line( size = 1) +
                              geom_ribbon(aes(ymin = lower__, ymax = upper__), color = NA, alpha = 0.1) +
                              theme_bw() +
                              labs(x = "Day", y="MDS1", title = "Day")
                            
                            plot_material <-
                              x$`material3` %>%
                              ggplot(aes(x = material3, y = estimate__, color = material3)) +
                              geom_point() +
                              geom_errorbar(aes(ymax = upper__, ymin = lower__), width = .4) +
                              theme_bw() +
                              labs(x = "Material", y="MDS1", color = "Material", fill = "Material", title = "Material") +
                              scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                            "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                            "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                            "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                            "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                            "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
                              
                            plot_interaction <-
                              x$`day:material3` %>%
                              ggplot(aes(x = day, y = estimate__, color = material3)) +
                              geom_line( size = 1) +
                              geom_ribbon(aes(ymin = lower__, ymax = upper__,  fill = material3), color = NA, alpha = 0.1) +
                              theme_bw() +
                              labs(x = "Day", y="MDS1", color = "Material", fill = "Material", title = "Interaction") +
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
                              
                            
                            #res <- plot_day/plot_material/plot_interaction
                            res <- plot_interaction
                            return(res)
                            
                          })) %>%
  mutate(plot2_mds2.1 = map(plot2_mds2,
                            function(x){
                              plot_day <-
                                x$`day` %>%
                                ggplot(aes(x = day, y = estimate__)) +
                                geom_line( size = 1) +
                                geom_ribbon(aes(ymin = lower__, ymax = upper__), color = NA, alpha = 0.1) +
                                theme_bw() +
                                labs(x = "Day", y="MDS1", title = "Day")
                              
                              plot_material <-
                                x$`material3` %>%
                                ggplot(aes(x = material3, y = estimate__, color = material3)) +
                                geom_point() +
                                geom_errorbar(aes(ymax = upper__, ymin = lower__), width = .4) +
                                theme_bw() +
                                labs(x = "Material", y="MDS1", color = "Material", fill = "Material", title = "Material") +
                                scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                              "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                              "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                              "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                              "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                              "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))
                              
                              plot_interaction <-
                                x$`day:material3` %>%
                                ggplot(aes(x = day, y = estimate__, color = material3)) +
                                geom_line( size = 1) +
                                geom_ribbon(aes(ymin = lower__, ymax = upper__,  fill = material3), color = NA, alpha = 0.1) +
                                theme_bw() +
                                labs(x = "Day", y="MDS1", color = "Material", fill = "Material", title = "Interaction") +
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
                              
                              #res <- plot_day/plot_material/plot_interaction
                              res <- plot_interaction
                              return(res)
                              
                            }))  %>%
  mutate(plot_mds1_day = map(df.nmds,
                             function(x){
                               x %>%
                                 ggplot(aes(x = day, y = MDS1, color = material3, group = id2)) +
                                 geom_point(alpha = 0.6) +
                                 geom_smooth(method = "lm", se =F, alpha = 0.6) +
                                 theme_bw() +
                                 labs(color = "Material", x = "Incubation Day") +
                                 scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                               "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                               "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                               "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                               "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                               "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
                             })) %>%
  mutate(plot_mds1_day.2 = map2(plot2_mds1.1, plot_mds1_day,
                                function(x, y){
                                  df.tmp1 <- x$data
                                  df.tmp2 <- y$data
                                  ggplot() +
                                    layer(data = df.tmp2,
                                          mapping = aes(x = day, y = MDS1, color = material3, group = id2, shape = test, fill = material3),
                                          geom = "point", stat = "identity", position = "identity",
                                          params = list(alpha = 0.6)) +
                                    layer(data = df.tmp1,
                                          mapping = aes(x = day, y = estimate__, color = material3),
                                          geom = "line", stat = "identity", position = "identity",
                                          params = list(size = 1)) +
                                    layer(data = df.tmp1,
                                          mapping = aes(x = day, y = estimate__, color = material3, ymin = lower__, ymax = upper__,  fill = material3), 
                                          geom = "ribbon", stat = "identity", position = "identity",
                                          params = list(color = NA, alpha = 0.1)) +
                                    theme_bw() +
                                    labs(x = "Day", y="nMDS1", color = "Material", fill = "Material") +
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
                                  
                                }))


  
for (i in 1:nrow(DF_PICRUSt_2)){
  db <- DF_PICRUSt_2[i,] %>% dplyr::select(db) %>% .[[1]]
  plot <- DF_PICRUSt_2[i,] %>% dplyr::select(plot.nmds) %>% .[[1]] %>% .[[1]]
  plot2 <- DF_PICRUSt_2[i,] %>% dplyr::select(plot_mds1_day) %>% .[[1]] %>% .[[1]]
  plot_mcmc_mds1 <- DF_PICRUSt_2[i,] %>% dplyr::select(mcmc_mds1) %>% .[[1]] %>% .[[1]]
  plot_mcmc_mds2 <- DF_PICRUSt_2[i,] %>% dplyr::select(mcmc_mds2) %>% .[[1]] %>% .[[1]]
  plot_marginal.effects_mds1 <- DF_PICRUSt_2[i,] %>% dplyr::select(plot2_mds1.1) %>% .[[1]] %>% .[[1]]
  plot_marginal.effects_mds2 <- DF_PICRUSt_2[i,] %>% dplyr::select(plot2_mds2.1) %>% .[[1]] %>% .[[1]]
  plot_marginal.effects_mds1.2 <- DF_PICRUSt_2[i,] %>% dplyr::select(plot_mds1_day.2) %>% .[[1]] %>% .[[1]]
  
  ggsave(plot = plot, file = str_c(CWD, "/res/picrust2/nMDS_",db, ".png", sep = ""), height = 4, width = 7.2)
  ggsave(plot = plot2, file = str_c(CWD, "/res/picrust2/MDS1vsDay_",db, ".png", sep = ""), height = 3, width = 5)
  ggsave(plot = plot_mcmc_mds1, file = str_c(CWD, "/res/picrust2/plot_mcmc_mds1_",db, ".png", sep = ""), height = 5, width = 5)
  ggsave(plot = plot_mcmc_mds2, file = str_c(CWD, "/res/picrust2/plot_mcmc_mds2_",db, ".png", sep = ""), height = 5, width = 5)
  ggsave(plot = plot_marginal.effects_mds1, file = str_c(CWD, "/res/picrust2/plot_marginal.effects_mds1_",db, ".png", sep = ""), height = 3, width = 5)
  ggsave(plot = plot_marginal.effects_mds2, file = str_c(CWD, "/res/picrust2/plot_marginal.effects_mds2_",db, ".png", sep = ""), height = 3, width = 5)
  ggsave(plot = plot_marginal.effects_mds1.2, file = str_c(CWD, "/res/picrust2/MDS1vsDay_bhm_", db, ".png", sep = ""), height = 3, width = 5)
  
  trace_mds1 <- DF_PICRUSt_2[i,] %>% dplyr::select(plot1_mds1) %>% .[[1]] %>% .[[1]]
  trace_mds2 <- DF_PICRUSt_2[i,] %>% dplyr::select(plot1_mds2) %>% .[[1]] %>% .[[1]]
  
  for(ii in 1:length(trace_mds1)){
    plot.tmp1 <- trace_mds1[[ii]]
    plot.tmp2 <- trace_mds2[[ii]]
    ggsave(plot = plot.tmp1, file = str_c(CWD, "/res/picrust2/trace_mds1_", ii,"_", db, ".png", sep = ""), height = 5, width = 5)
    ggsave(plot = plot.tmp2, file = str_c(CWD, "/res/picrust2/trace_mds2_", ii,"_", db, ".png", sep = ""), height = 5, width = 5)
  }
  
}


#write_rds(DF_PICRUSt_2, str_c(Dir_res_integral, "/PICRUSt2_output.rds"))
write_rds(DF_PICRUSt_2, str_c(Dir_res_integral, "/PICRUSt2_output.rds"), compress = "gz")

Dir_res_picrust2_g1 <- str_c(CWD, "/res/picrust2/sp.score.g1")
Dir_res_picrust2_g2 <- str_c(CWD, "/res/picrust2/sp.score.g2")
Dir_res_picrust2_g3 <- str_c(CWD, "/res/picrust2/sp.score.g3")

if(!file.exists(Dir_res_picrust2_g1)){
  dir.create(Dir_res_picrust2_g1)
}
if(!file.exists(Dir_res_picrust2_g2)){
  dir.create(Dir_res_picrust2_g2)
}
if(!file.exists(Dir_res_picrust2_g3)){
  dir.create(Dir_res_picrust2_g3)
}


KO.SCORE <-
  DF_PICRUSt_2 %>%
  filter(db == "data_ko.2") %>%
  dplyr::select(nmds) %>%
  .[[1]] %>% 
  .[[1]] %>% 
  .$species %>%
  as_tibble(rownames = "KO")

Group.1_list <- unique(KO_list$Group.1)
Group.2_list <- unique(KO_list$Group.2)
Group.3_list <- unique(KO_list$Group.3)

all.score.kegg <- as_tibble()

for(i in 1:length(Group.3_list)){
  tmp <-
    KO_list %>%
    filter(Group.3 == Group.3_list[i])  %>%
    left_join(KO.SCORE, ., by = "KO") 
  
  tmp2 <-
    tmp %>%
    filter(!is.na(Group.3))
  
  all.score.kegg <-
    bind_rows(all.score.kegg, tmp2)
  
  plot_org <-
    ggplot() +
    geom_point(data = tmp, aes(x = MDS1, y = MDS2), color = "grey") +
    theme_bw() +
    ggtitle(Group.3_list[i]) 
  
  plot1 <-
    plot_org +
    geom_point(data = tmp %>% filter(!is.na(Group.3)), aes(x = MDS1, y = MDS2), color = "#F8766D") 
  
  if(tmp %>% filter(!is.na(Group.3)) %>% nrow() >2){
    plot2 <-
      plot_org +
      stat_density_2d(data = tmp %>% filter(!is.na(Group.3)), geom = 'polygon', aes(x = MDS1, y =MDS2, alpha = ..level.., fill = Group.3)) +
      theme(legend.position = "none")
  }  else {
    plot2 <- plot1
  }
  
  plot3 <- 
    bind_rows(tmp %>% mutate(Group.3 = "all"), tmp %>% filter(!is.na(Group.3)))  %>%
    mutate(Group.3 = factor(Group.3, levels = c("all", Group.3_list[i]))) %>%
    ggplot(aes(x = MDS1, fill = Group.3)) +
    geom_density(alpha = 0.4) + theme_bw() + scale_fill_manual(values = c("grey", "#00BFC4")) +
    ggtitle(Group.3_list[i]) +
    theme(legend.position = "none")
    

  ggsave(plot = plot1,
         file = str_c(Dir_res_picrust2_g3, "/sp.score.g3_", i, ".png", sep = ""),
         height = 5, width = 5)
  ggsave(plot = plot2,
         file = str_c(Dir_res_picrust2_g3, "/density_sp.score.g3_", i, ".png", sep = ""),
         height = 5, width = 5)
  ggsave(plot = plot3,
         file = str_c(Dir_res_picrust2_g3, "/density_mds1.g3_", i, ".png", sep = ""),
         height = 2, width = 4)
  
}

all.score.kegg.2 <-
  all.score.kegg %>%
  nest(-Group.3) %>%
  mutate(median.mds1 = map_dbl(data, ~median(.$MDS1)),
         num.ko = map_dbl(data, ~nrow(.))) %>%
  arrange(median.mds1) %>%
  mutate(Group.3 = factor(Group.3, levels = .$Group.3)) %>%
  filter(num.ko > 40) %>%
  mutate(rank = row_number(), rank2 = nrow(.) - rank + 1) %>%
  filter(rank <= 10 | rank2 <= 10) %>%
  unnest()  %>%
  nest() %>%
  mutate(multcomp = map(data, 
                        function(x){
                          
                          glht(aov(MDS1 ~ Group.3, data = x), linfct = mcp(Group.3 = "Tukey")) %>% 
                            summary() %>% 
                            cld(decreasing = T) %>%
                            .$mcletters %>% .$Letters %>% as_tibble(rownames = "Group.3")
                        })) %>%
  mutate(plot = map2(data, multcomp,
                     function(x, y){
                       ggplot() +
                         geom_boxplot(data = x, aes(x = Group.3, y = MDS1)) +
                         geom_text(data = y, aes(x = Group.3, label = value), y = 1.5) + 
                         coord_flip() + 
                         ylim(-1.2, 1.5) + 
                         theme_bw() + 
                         xlab("KEGG Orthology")
                     })) %>%
  mutate(multcomp.2 = map(data,
                                  function(x){
                                    Group.3.order <-  x$Group.3 %>% unique() %>% sort()
                                    
                                    tmp1 <-
                                      x %>%
                                      nest() %>%
                                      mutate(multcomp = map(data,
                                                            function(x){
                                                              #tic()
                                                              res <- my_own_nparcomp(MDS1 ~ Group.3, data=x, asy.method = "mult.t",
                                                                                     type = "Tukey",alternative = "two.sided",info = FALSE)
                                                              #toc()
                                                              return(res)
                                                            })) %>%
                                      mutate(letters = map(multcomp, func_letter_nparcomp))
                                    
                                    
                                    data.letters <- 
                                      tmp1 %>%
                                      dplyr::select( letters) %>%
                                      unnest() %>%
                                      rename(Group.3 = "param") %>%
                                      mutate(Group.3 = factor(Group.3, levels = Group.3.order))
                                    
                                    data.out <-
                                      tmp1 %>%
                                      dplyr::select( data) %>%
                                      unnest() %>%
                                      mutate(Group.3 = factor(Group.3, levels = Group.3.order))
                                    
                                    list(data.out, data.letters)  
                                  })) %>%
  mutate(plot.2 = map(multcomp.2,
                      function(x){
                        data <- x[[1]] %>% mutate(position = ifelse(rank <= 10, "lower", ifelse(rank2 <= 10, "higher", "other"))) 
                        df.letters <- x[[2]]
                        
                        ggplot() +
                          geom_jitter(data = data, aes(x = Group.3, y = MDS1, color = position), alpha = 0.2) +
                          geom_boxplot(data = data, aes(x = Group.3, y = MDS1), fill = NA, outlier.colour = NA) +
                          geom_text(data = df.letters, aes(x = Group.3, y = 1.15, label = fin.letter)) +
                          theme_bw() +  
                          
                          coord_flip() + 
                          ylim(-1.2, 1.5) + 
                          labs(x = "Pathway", y = "MDS1", color = "Position on MDS1") #+
                          #theme(legend.position = "bottom")
                        
                      })) 


ggsave(plot = all.score.kegg.2 %>% dplyr::select(plot.2) %>% .[[1,1]] %>% .[[1]],
       file = str_c(Dir_res_picrust2, "/KO_on_nmds1.pdf"),
       width = 8, height = 8)

select.list <-
  c(Group.3_list[str_detect(Group.3_list, "Biofilm")], 
    Group.3_list[str_detect(Group.3_list, "Cell growth")],
    Group.3_list[str_detect(Group.3_list, "Flagellar")])
  
  c(Group.3_list[str_detect(Group.3_list, "Biofilm")], 
    Group.3_list[str_detect(Group.3_list, "Cell growth")])




DF_Density <-as_tibble()
for(i in 1:length(select.list)){
  tmp <-
    KO_list %>%
    filter(Group.3 == select.list[i])  %>%
    left_join(., KO.SCORE, by = "KO") %>%
    mutate(Group.4 = ifelse(str_detect(Group.3, "Biofilm"), "Biofilm formation",
                            ifelse(str_detect(Group.3, "Cell"), "Cell growth", 
                                   ifelse(str_detect(Group.3, "Flagellar"), "Flagellar assembly", "Other"))))
  DF_Density <-
    bind_rows(DF_Density, tmp) 
}

DF_Density <-
  KO_list %>%
  left_join(., KO.SCORE, by = "KO") %>%
  mutate(Group.4 = "all") %>%
  bind_rows(., DF_Density) 

dplot_axis1 <-
  DF_Density %>%
  ggplot(aes(x = MDS1, fill = Group.4)) +
  geom_density(alpha = 0.4) +
  theme_bw() +
  scale_fill_manual(values = c("grey", hue_pal()(3)[1], hue_pal()(3)[2],  hue_pal()(3)[3])) +
  theme(legend.position = "none",
        axis.title= element_blank(),
        axis.text = element_blank()) 

nmds_select <-
  ggplot() +
  geom_point(data = DF_Density %>% filter(Group.4 == "all") %>% 
               mutate(Gene.group = ifelse(str_detect(Group.3, "Biofilm"), "Biofilm formation", 
                                          ifelse(str_detect(Group.3, "Cell growth"), "Cell growth",
                                                 ifelse(str_detect(Group.3, "Flagellar assembly"), "Flagellar assembly",  "All")))),
             aes(x = MDS1, y = MDS2, color = Gene.group)) +   scale_color_manual(values = c("grey",  hue_pal()(3)[1], hue_pal()(3)[2],  hue_pal()(3)[3])) +
  geom_point(data = DF_Density %>% filter(Group.4 == "Biofilm formation"),
             aes(x = MDS1, y = MDS2), color = hue_pal()(3)[1]) +
  geom_point(data = DF_Density %>% filter(Group.4 == "Cell growth"),
             aes(x = MDS1, y = MDS2), color = hue_pal()(3)[2]) +
  geom_point(data = DF_Density %>% filter(Group.4 == "Flagellar assembly"),
             aes(x = MDS1, y = MDS2), color = hue_pal()(3)[3]) +
  theme_bw() +
  theme(legend.position = c(0.17, 0.84),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.text = element_text(colour = "black", size = 8)
        ) +
  labs(x = "nMDS1", y = "nMDS2")

library(patchwork)
plot_nmds_select <- 
  dplot_axis1 + 
  nmds_select + 
  plot_layout(ncol = 1, heights = c(1, 4))

ggsave(plot = plot_nmds_select, 
       file = str_c(CWD, "/res/picrust2/nMDS_biofilm_growth.png", sep = ""), height = 6, width = 5)

# # sumamrize KO into Pathway for bn analysis
# 
# DF_each_KO <-
#   DF_PICRUSt_2 %>%
#   filter(db == "data_ko.2") %>%
#   dplyr::select(data3) %>%
#   unnest() %>%
#   gather(-c(test, material, material2, material3, time, bs, day), key = "KO", value = "prop") %>%
#   nest(-KO) %>%
#   left_join(., KO_list %>% dplyr::select(KO, gene.name, description) %>% distinct(), by = "KO") %>%
#   mutate(anova_4way = map(data,
#                           function(x){
#                             summary(aov(prop * 100 ~ test * material3 * day * bs, data = x))
#                           })) %>%
#   mutate(anova_4way_res = map_chr(anova_4way,
#                                   function(x){
#                                     tmp1 <- x[[1]] %>% rownames() %>% as_tibble()
#                                     tmp2 <- x[[1]] %>%  as_tibble()
#                                     tmp3 <-
#                                       bind_cols(tmp1, tmp2) %>%
#                                       set_names(c("factor", "df", "Sum_seq", "mean_seq", "F.value", "p.value")) %>%
#                                       mutate(factor = str_replace(factor, pattern = "test", replacement = "Test"),
#                                              factor = str_replace(factor, pattern = "material3", replacement = "Material"),
#                                              factor = str_replace(factor, pattern = "day", replacement = "Day"),
#                                              factor = str_replace(factor, pattern = "bs", replacement = "Biofilm"),
#                                              factor = str_trim(factor))
#                                     tmp4 <- 
#                                       tmp3 %>%
#                                       mutate(factor.num = str_split(factor, pattern = ":"),
#                                              factor.num = map_int(factor.num, ~length(.))) %>%
#                                       nest(-factor.num) %>%
#                                       mutate(tag = map_chr(data,
#                                                            function(x){
#                                                              x %>%
#                                                                filter(p.value < 0.05) %>%
#                                                                mutate(factor = str_c(factor, "*")) %>%
#                                                                .$factor %>% 
#                                                                str_c(collapse = ",") %>%
#                                                                str_c("    ", ., sep = "")
#                                                            })) 
#                                     
#                                     res <- tmp4$tag %>% str_c(collapse = ("\n"))
#                                     
#                                     return(res)
#                                   })) 
#   #filter(str_detect(anova_4way_res, "Material") & !str_detect(anova_4way_res, ":Material") & !str_detect(anova_4way_res, "Material:")) %>%
#   mutate(plot = pmap(list(data, anova_4way_res), 
#                    function(XX, YY){
#                      XX %>%
#                        ggplot(aes(x = material2, y = prop * 100, color = material3)) + 
#                        geom_boxplot() + 
#                        theme_bw() + 
#                        ylab("Proportion (%)") + 
#                        xlab("Material") +
#                        labs(color = "Material") +
#                        ggtitle(str_c("4-way ANOVA;\n", YY)) +
#                        facet_grid(test ~ .)  +
#                        scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                      "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                      "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                      "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                      "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                      "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) 
#                    })) %>%
#   mutate(plot2 = pmap(list(data, KO, gene.name, description), 
#                       function(XX, x, y, z){
#                         XX %>%
#                           mutate(id = str_c(test, material, bs, sep = "_")) %>%
#                           mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
#                           ggplot(aes(x = day, y = 100* prop, color = material3, group = id, shape = bs2)) +
#                           geom_line(alpha = 0.4) +
#                           geom_point(alpha = 0.4) +
#                           #geom_smooth(se = F, alpha = 0.4) +
#                           facet_grid(test ~ .) +
#                           theme_bw() +
#                           ylab("Proportion (%)") + 
#                           xlab("Incubation Day") +
#                           labs(color = "Material", shape = "Biofilm/Surface") +
#                           ggtitle(str_c(x, " ", y))+
#                           scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                         "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                         "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                         "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                         "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                         "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]))
#                       })) %>%
#   mutate(plot3 = pmap(list(plot, plot2, KO, gene.name, description, 
#                            anova_4way_res), 
#                       function(XX, YY, x, y, z, a){
#                         XX + theme(legend.position = "none") + theme(plot.title = element_blank()) + 
#                           YY + theme(plot.title = element_blank()) +
#                           plot_annotation(title = x,
#                                           subtitle = str_c(y,"\n", z,"\n",
#                                                            "4-way ANOVA;\n", a))
#                       })) %>%
#   mutate(t.test = map(data, 
#                       function(x){
#                         t.test(x %>% filter(material2 == "PHBH") %>% .$prop, 
#                                x %>% filter(material2 != "PHBH") %>% .$prop, 
#                                var.equal = T)})) %>%
#   mutate(t.test.p = map_dbl(t.test, function(x){x$p.value})) %>%
#   arrange(t.test.p)
#   
# 
# Dir_res_picrust2_ts<- str_c(CWD, "/res/picrust2/ts")
# 
# if(!file.exists(Dir_res_picrust2_ts)){
#   dir.create(Dir_res_picrust2_ts)
# }
# 
# 
# for (i in 1:nrow(DF_each_KO)){
#   ko.tmp <-  DF_each_KO[i,] %>% dplyr::select(KO) %>% .[[1]]
#   gene.tmp <- DF_each_KO[i,] %>% dplyr::select(gene.name) %>% .[[1]]
#   desc.tmp <-  DF_each_KO[i,] %>% dplyr::select(description) %>% .[[1]]
#   plot <- DF_each_KO[i,] %>% dplyr::select(plot3) %>% .[[1]] %>% .[[1]]
#   
#   ggsave(plot =  plot, 
#          file = str_c(Dir_res_picrust2_ts, "/combine_plot_", i ,"_", ko.tmp, ".png"), height = 10, width = 10)
# }


#### for GLMM
# K05973: phaZ
# K08095: cutinase
# K01046: triacylglycerol lipase
# K12298: bile salt-stimulated lipase 
# K14674: TAG lipase / steryl ester hydrolase / phospholipase A2 / LPA acyltransferase
# K01188, K05349, K05350: b-glucosidase
KO_target <-
  DF_PICRUSt_2 %>%
  filter(db == "data_ko.2") %>%
  mutate(data3 = map(data3, function(x){x %>% dplyr::select(c(test, material, material2, material3, bs, day, K05973, K08095, K01046, K01188, K07406, K12309, K11935, K11937))})) %>%
  dplyr::select(data3) %>%
  unnest() %>%
  gather(-c(test, material, material2, material3, bs, day), key = "KO", value = "prop") %>%
  mutate(bs2 = ifelse(bs == "bio"|bs == "not", "Biofilm", ifelse(bs == "surf", "Surface", "others"))) %>%
  mutate(prop = 100 * prop,
         id2 = str_c(test, material, bs, sep = "_"),
         id3 = str_c(test, material, sep = "_")) %>%
  nest(-KO)

#write_rds(KO_target, "KO_target.rds")

# 2 factors
KO_target_2 <-
  KO_target %>%
  left_join(., KO_list %>% dplyr::select(KO, gene.name, description) %>% distinct(), by = "KO") %>%
  mutate(glmm = map(data,
                    function(x){
                      brm(formula = prop ~ day * material3 + (1|id2),
                          family = gaussian(link = "identity"),
                          data = x,
                          seed = 1,
                          iter = 2000,
                          warmup = 1000)
                    })) %>%
  mutate(plot1 = map(glmm, function(x){plot(x)})) %>%
  mutate(fixed_effect = map(glmm,
                            function(x){
                              x %>%
                                summary() %>%
                                .$fixed %>%
                                as_tibble(rownames = "factor")
                            })) %>%
  mutate(mcmc = pmap(list(glmm, KO, gene.name, description),
                     function(x, a, b, c){
                       x %>%
                         spread_draws(b_Intercept, b_day,`b_material3PHBH_HH6%`, b_material3PCL, b_material3PBSA, b_material3PBS, b_material3PBAT, 
                                      `b_day:material3PHBH_HH6%`,`b_day:material3PCL`, `b_day:material3PBSA`, `b_day:material3PBS`, `b_day:material3PBAT`) %>%
                         gather(-c(.chain, .iteration, .draw), key = "factor", value = "val") %>%
                         mutate(factor = str_replace(factor, pattern = "material3", replace = "material")) %>%
                         mutate(factor2 = ifelse(str_detect(factor, "b_day:material"), "Interaction", 
                                                 ifelse(str_detect(factor, "b_material"), "Material", 
                                                        ifelse(factor == "b_day", "Day",
                                                               ifelse(factor == "b_Intercept", "Intercept", "other"))))) %>%
                         mutate(factor2 = factor(factor2, levels = c("Intercept", "Day", "Material", "Interaction"))) %>%
                         mutate(factor = str_replace(factor, pattern = "b_", replacement = "")) %>% 
                         mutate(factor = str_replace(factor, pattern = "material", replacement = "")) %>%
                         mutate(factor = str_replace(factor, pattern = "day", replacement = "Day")) %>%
                         ggplot(aes(x = val, y = factor, fill = factor(stat(quantile)))) +
                         stat_density_ridges(geom = "density_ridges_gradient",
                                             calc_ecdf = TRUE,
                                             quantiles = c(0.025, 0.975)) +
                         scale_fill_manual(name = "Probability", 
                                           values = c("#A0A0A0A0", "#00BFC4", "#A0A0A0A0"),
                                           labels = c("0 - 2.5 %", "2.5 - 97.5 %", "97.5 - 100 %")) +
                         facet_grid2(factor2 ~ ., scales = "free", independent = "x") +
                         theme_classic() +
                         geom_vline(xintercept = 0, color = "#F8766D", linetype = 2, size = 1) +
                         theme( strip.background = element_blank()) +
                         xlab("Estimate") +
                         ylab("Factor") +
                         labs(title = a,
                              subtitle = str_c(b, "\n", c))
                     })) %>%
  mutate(plot2 = map(glmm,
                     function(x){
                       marginal_effects(x)
                     })) %>%
  mutate(plot2.1 = pmap(list(plot2, KO, gene.name, description),
                            function(x, a, b, c){
                              #col.pal <- hue_pal()(5)
                              plot_day <-
                                x$`day` %>%
                                ggplot(aes(x = day, y = estimate__)) +
                                geom_line( size = 1) +
                                geom_ribbon(aes(ymin = lower__, ymax = upper__), color = NA, alpha = 0.1) +
                                theme_bw() +
                                labs(x = "Day", y="Proportion", title = "Day")
                              
                              plot_material <-
                                x$`material3` %>%
                                ggplot(aes(x = material2, y = estimate__, color = material3)) +
                                geom_point() +
                                geom_errorbar(aes(ymax = upper__, ymin = lower__), width = .4) +
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
                                                             "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                labs(x = "Material", y="Proportion", color = "Material", title = "Material") 
                              
                              plot_interaction <-
                                x$`day:material3` %>%
                                ggplot(aes(x = day, y = estimate__, color = material3)) +
                                geom_line( size = 1) +
                                geom_ribbon(aes(ymin = lower__, ymax = upper__,  fill = material3), color = NA, alpha = 0.1) +
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
                                                             "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
                                labs(x = "Day", y="Proportion", color = "Material", fill = "Material", title = "Interaction", shape = "Biofilm/Surface")
                              
                              #res <- 
                              #    plot_day/plot_material/plot_interaction +
                              #    plot_annotation(title = a,
                              #                    subtitle = str_c(b, "\n", c))
                              res <- plot_interaction +
                                plot_annotation(title = a, subtitle = str_c(b, "\n", c))
                              
                              return(res)
                              
                            })) %>%  mutate(plot2.2 = pmap(list(data, plot2.1, KO, gene.name, description),
                                                           function(x, y, a, b, c){
                                                             df.tmp1 <- x
                                                             df.tmp2 <- y$data
                                                             ggplot() +
                                                               layer(data = df.tmp1,
                                                                     mapping = aes(x = day, y = prop, color = material3, group = id2, shape = bs2),
                                                                     geom = "point", stat = "identity", position = "identity",
                                                                     params = list(alpha = 0.6)) +
                                                               layer(data = df.tmp2,
                                                                     mapping = aes(x = day, y = estimate__, color = material3),
                                                                     geom = "line", stat = "identity", position = "identity",
                                                                     params = list(size = 1)) +
                                                               layer(data = df.tmp2,
                                                                     mapping = aes(x = day, y = estimate__, color = material3, ymin = lower__, ymax = upper__,  fill = material3), 
                                                                     geom = "ribbon", stat = "identity", position = "identity",
                                                                     params = list(color = NA, alpha = 0.1)) +
                                                               theme_bw() +
                                                               labs(x = "Day", y="Proportion (%)", color = "Material", fill = "Material", shape = "Biofilm/Surface") +
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
                                                               plot_annotation(title = a, subtitle = str_c(b, "\n", c))
                                                             
                                                           }))



# KO_target_2.1 <-
#   KO_target %>%
#   left_join(., KO_list %>% dplyr::select(KO, gene.name, description) %>% distinct(), by = "KO") %>%
#   mutate(data = map(data, function(x){x %>% filter(material3 != "PBAT" & material3 != "PBS"& material3 != "PHBH_HH6%")})) %>%
#   mutate(glmm = map(data,
#                     function(x){
#                       brm(formula = prop ~ day * material3 + (1|id2),
#                           family = gaussian(link = "identity"),
#                           data = x,
#                           seed = 1,
#                           iter = 3000,
#                           warmup = 2000)
#                     })) %>%
#   mutate(plot1 = map(glmm, function(x){plot(x)})) %>%
#   mutate(fixed_effect = map(glmm,
#                             function(x){
#                               x %>%
#                                 summary() %>%
#                                 .$fixed %>%
#                                 as_tibble(rownames = "factor")
#                             })) %>%
#   mutate(mcmc = pmap(list(glmm, KO, gene.name, description),
#                      function(x, a, b, c){
#                        x %>%
#                          spread_draws(b_Intercept, b_day, b_material3PCL, b_material3PBSA, `b_day:material3PCL`, `b_day:material3PBSA`) %>%
#                          gather(-c(.chain, .iteration, .draw), key = "factor", value = "val") %>%
#                          mutate(factor = str_replace(factor, pattern = "material3", replace = "material")) %>%
#                          mutate(factor2 = ifelse(str_detect(factor, "b_day:material"), "Interaction", 
#                                                  ifelse(str_detect(factor, "b_material"), "Material", 
#                                                         ifelse(factor == "b_day", "Day",
#                                                                ifelse(factor == "b_Intercept", "Intercept", "other"))))) %>%
#                          mutate(factor2 = factor(factor2, levels = c("Intercept", "Day", "Material", "Interaction"))) %>%
#                          mutate(factor = str_replace(factor, pattern = "b_", replacement = "")) %>% 
#                          mutate(factor = str_replace(factor, pattern = "material", replacement = "")) %>%
#                          mutate(factor = str_replace(factor, pattern = "day", replacement = "Day")) %>%
#                          ggplot(aes(x = val, y = factor, fill = factor(stat(quantile)))) +
#                          stat_density_ridges(geom = "density_ridges_gradient",
#                                              calc_ecdf = TRUE,
#                                              quantiles = c(0.025, 0.975)) +
#                          scale_fill_manual(name = "Probability", 
#                                            values = c("#A0A0A0A0", "#00BFC4", "#A0A0A0A0"),
#                                            labels = c("0 - 2.5 %", "2.5 - 97.5 %", "97.5 - 100 %")) +
#                          facet_grid2(factor2 ~ ., scales = "free", independent = "x") +
#                          theme_classic() +
#                          geom_vline(xintercept = 0, color = "#F8766D", linetype = 2, size = 1) +
#                          theme( strip.background = element_blank()) +
#                          xlab("Estimate") +
#                          ylab("Factor") +
#                          labs(title = a,
#                               subtitle = str_c(b, "\n", c))
#                      })) %>%
#   mutate(plot2 = map(glmm,
#                      function(x){
#                        marginal_effects(x)
#                      })) %>%
#   mutate(plot2.1 = pmap(list(plot2, KO, gene.name, description),
#                        function(x, a, b, c){
#                          col.pal <- hue_pal()(5)
#                          plot_day <-
#                            x$`day` %>%
#                            ggplot(aes(x = day, y = estimate__)) +
#                            geom_line( size = 1) +
#                            geom_ribbon(aes(ymin = lower__, ymax = upper__), color = NA, alpha = 0.1) +
#                            theme_bw() +
#                            labs(x = "Day", y="Proportion (%)", title = "Day")
#                          
#                          plot_material <-
#                            x$`material3` %>%
#                            ggplot(aes(x = material3, y = estimate__, color = material3)) +
#                            geom_point() +
#                            geom_errorbar(aes(ymax = upper__, ymin = lower__), width = .4) +
#                            theme_bw() +
#                            scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                          "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                          "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                          "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                          "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                          "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                            labs(x = "Material", y="Proportion (%)", color = "Material", title = "Material") 
#                          
#                          plot_interaction <-
#                            x$`day:material3` %>%
#                            ggplot(aes(x = day, y = estimate__, color = material3)) +
#                            geom_line( size = 1) +
#                            geom_ribbon(aes(ymin = lower__, ymax = upper__,  fill = material3), color = NA, alpha = 0.1) +
#                            theme_bw() +
#                            scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                          "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                          "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                          "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                          "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                          "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                            scale_fill_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                         "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                         "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                         "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                         "PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                         "PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16])) +
#                            labs(x = "Day", y="Proportion (%)", color = "Material", fill = "Material", title = "Interaction")
#                          
#                          #res <- 
#                          #    plot_day/plot_material/plot_interaction +
#                          #    plot_annotation(title = a,
#                          #                    subtitle = str_c(b, "\n", c))
#                          res <- plot_interaction +
#                            plot_annotation(title = a, subtitle = str_c(b, "\n", c))
#                          
#                          return(res)
#                          
#                        })) %>%
#   mutate(plot2.2 = pmap(list(data, plot2.1, KO, gene.name, description),
#                         function(x, y, a, b, c){
#                           df.tmp1 <- x
#                           df.tmp2 <- y$data
#                           ggplot() +
#                             layer(data = df.tmp1,
#                                   mapping = aes(x = day, y = prop, color = material3, group = id2, shape = bs2),
#                                   geom = "point", stat = "identity", position = "identity",
#                                   params = list(alpha = 0.6)) +
#                             layer(data = df.tmp2,
#                                   mapping = aes(x = day, y = estimate__, color = material3),
#                                   geom = "line", stat = "identity", position = "identity",
#                                   params = list(size = 1)) +
#                             layer(data = df.tmp2,
#                                   mapping = aes(x = day, y = estimate__, color = material3, ymin = lower__, ymax = upper__,  fill = material3), 
#                                   geom = "ribbon", stat = "identity", position = "identity",
#                                   params = list(color = NA, alpha = 0.1)) +
#                             theme_bw() +
#                             labs(x = "Day", y="Proportion (%)", color = "Material", fill = "Material", shape = "Biofilm/Surface") +
#                             scale_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                           #"PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                           "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                           "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                           #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                           #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]
#                                                           )) +
#                             scale_fill_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
#                                                          #"PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
#                                                          "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
#                                                          "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
#                                                          #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
#                                                          #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16]
#                                                          )) +
#                             plot_annotation(title = a, subtitle = str_c(b, "\n", c))
#                           
#                         }))



for(i in 1:length(KO_target_2)){
  KO.tmp <- KO_target_2[i,] %>% dplyr::select(KO) %>% .[[1]]
  plot.tmp1 <- KO_target_2[i,] %>% dplyr::select(mcmc) %>% .[[1]] %>% .[[1]]
  plot.tmp2 <- KO_target_2[i,] %>% dplyr::select(plot2.1) %>% .[[1]] %>% .[[1]]
  plot.tmp3 <- KO_target_2[i,] %>% dplyr::select(plot2.2) %>% .[[1]] %>% .[[1]]
  # plot.tmp1.1 <- KO_target_2.1[i,] %>% select(mcmc) %>% .[[1]] %>% .[[1]]
  # plot.tmp2.1 <- KO_target_2.1[i,] %>% select(plot2.1) %>% .[[1]] %>% .[[1]]
  # plot.tmp3.1 <- KO_target_2.1[i,] %>% select(plot2.2) %>% .[[1]] %>% .[[1]]
  
  ggsave(plot = plot.tmp1, file = str_c(CWD, "/res/picrust2/plot_mcmc_", i, "_", KO.tmp, ".png", sep = ""), height = 5, width = 5)
  ggsave(plot = plot.tmp2, file = str_c(CWD, "/res/picrust2/plot_marginal.effects_", i, "_",  KO.tmp, ".png", sep = ""), height = 4.5, width = 6)
  ggsave(plot = plot.tmp3, file = str_c(CWD, "/res/picrust2/plot_marginal.effects.1.2_", i, "_",  KO.tmp, ".png", sep = ""), height = 4.5, width = 6 )
  
  # ggsave(plot = plot.tmp1.1, file = str_c(CWD, "/res/picrust2/plot_mcmc.1_", i, "_", KO.tmp, ".png", sep = ""), height = 5, width = 5)
  # ggsave(plot = plot.tmp2.1, file = str_c(CWD, "/res/picrust2/plot_marginal.effects.1_", i, "_",  KO.tmp, ".png", sep = ""), height = 4.5, width = 6)
  # ggsave(plot = plot.tmp3.1, file = str_c(CWD, "/res/picrust2/plot_marginal.effects.1.2_", i, "_",  KO.tmp, ".png", sep = ""), height = 4.5, width = 6)
  
}
write_rds(KO_target_2, str_c(Dir_res_picrust2, "/KO_target_2factor_mcmc.rds"))
#write_rds(KO_target_2.1, str_c(Dir_res_picrust2, "/KO_target_2.1factor_mcmc.rds"))

######################
### Figs for thesis
######################

Figure.4_1 <- DF_PICRUSt_2[[1,5]][[1]] + theme(legend.position = "none")
Figure.4_2 <- DF_PICRUSt_2[[1,20]][[1]]+ theme(legend.position = "none") +   scale_x_continuous(breaks = seq(0,60,10)) 
legend <- DF_PICRUSt_2[[1,5]][[1]] %>% get_legend() %>% as_ggplot()
Figure.4_3 <- all.score.kegg.2 %>% dplyr::select(plot.2) %>% .[[1]] %>% .[[1]] + theme(axis.text.y=element_text(colour = c(rep(hue_pal()(2)[2],10), rep(hue_pal()(2)[1], 10))))
Figure.4 <-
  gridExtra::grid.arrange(Figure.4_1 + labs(title = "(A)") + theme(plot.title.position = "plot"), 
                          Figure.4_2 + labs(title = "(B)") + theme(plot.title.position = "plot"),
                          legend, 
                          Figure.4_3+ labs(title = "(C)") + theme(plot.title.position = "plot"), 
                          layout_matrix = rbind(c(1,2,3), c(4,4,4)))
ggsave(plot = Figure.4,
       file = str_c(Dir_res_final_fig, "/Figure.4.pdf", sep = ""), 
       height = 6, width = 10)

