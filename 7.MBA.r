rm(list = ls())

library(tidyverse)
library(arules) 
library(arulesViz)
library(tidygraph)
library(ggraph)
library(scales)
library(ggnewscale)
library(colorspace)

CWD <- getwd()
DIR_Integrate <- str_c(CWD, "/res/integral")
Dir_res_final_fig <- str_c(CWD, "/res/final.fig")
if(!file.exists(Dir_res_final_fig)){
  dir.create(Dir_res_final_fig)
}

DF_mba_mic <- read_rds(str_c(DIR_Integrate, "/DF_mba_mic.rds"))
#DF_mba_ko <- read_rds(str_c(DIR_Integrate, "/DF_mba_ko.rds"))
DF_combine <- read_rds(str_c(DIR_Integrate, "/DF_combine.rds"))

### 3 categories (experimentaal.factor, microbiome, metabolome)
MBA <-
  DF_combine %>%
  filter(param == "mic") %>% #omit picrust
  mutate(data2 = map(data, # top25%を1
                     function(x){
                       df2.1 <-
                         x %>%
                         dplyr::select(test, material3, day, bs) %>%
                         mutate(id = row_number()) %>% 
                         mutate(day = ifelse(day <= 4, "day.0-4", ifelse(day <= 14, "day.4-14", "day.14-"))) %>%
                         gather(-id, key = "exist", value = "param") %>% 
                         mutate(exist = T) %>% 
                         spread(key =param, value = exist) %>% 
                         arrange(id) %>%
                         mutate_all(~ifelse(is.na(.), F, .)) 
                       
                       df2.2 <-
                         x %>%
                         mutate(id = row_number()) %>%
                         dplyr::select(-c(test, material2, material3, day, bs)) %>%
                         select_if(function(x){length(unique(x)) > 1}) %>%
                         gather(-id, key = "param", value = "val") %>%
                         nest(-param) %>%
                         mutate(summary = map(data, function(x){summary(x$val)})) %>%
                         mutate(data2 = map2(data, summary, function(x,y){x %>% mutate(val2 = ifelse(val > y[[5]], T,  F))})) %>%
                         dplyr::select(-data, -summary) %>%
                         unnest(data2) %>%
                         dplyr::select(-val) %>%
                         spread(key = param, value = val2)
                       
                       DF_market.basket <-
                         left_join(df2.1, df2.2, by = "id") %>%
                         dplyr::select(-id) %>%
                         select_if(function(x){length(unique(x)) > 1}) 
                       
                       return(DF_market.basket)
                       
                     })) %>%
  mutate(transaction = map(data2, function(x){as(x, "transactions")})) %>%
  mutate(rules = map(transaction, function(x){apriori(x, parameter=list(support=0.063,confidence=0.25,maxlen=2))})) %>%
  mutate(rules_lis = map(rules, function(x){capture.output(inspect(x))})) %>%
  mutate(rules_table = map(rules_lis,
                           function(x){
                             res <-
                               x[-1]  %>%
                               as_tibble() %>%
                               mutate(value = str_replace_all(value, "; ", ";")) %>%
                               mutate(tmp = str_split(value, pattern = " ")) %>%
                               mutate(col.list = map(tmp, function(x){x[x != ""] %>% t()}))  %>% 
                               dplyr::select(col.list) %>%
                               unnest(col.list) %>%
                               .[[1]] %>%
                               as_tibble() %>%
                               set_names(c("id", "source", "direction", "target", "support", "confidence", "converge", "lift", "count")) %>%
                               dplyr::select(-direction, -id) %>%
                               type_convert() %>%
                               mutate(G.source = ifelse(str_detect(source, "d__"), "Microbiome",
                                                        ifelse(str_detect(source, "Metabolite"), "Metabolome", "Experimental.Factor")),
                                      G.target = ifelse(str_detect(target, "d__"), "Microbiome",
                                                        ifelse(str_detect(target, "Metabolite"), "Metabolome", "Experimental.Factor"))) %>% filter(!(source == "{}" | target == "{}")) %>%
                               mutate(link.mic.met = ifelse((G.source == "Metabolome" & G.target == "Microbiome")|(G.source == "Microbiome" & G.target == "Metabolome"), T, F)) %>%
                               mutate(link.with.ex = ifelse(G.source == "Experimental.Factor"|G.target == "Experimental.Factor", T, F)) %>%
                               mutate(link.with.PHBH10 = ifelse(source == "{PHBH_HH10%}"|target == "{PHBH_HH10%}", T, F)) %>%
                               mutate(link.with.PHBH6 = ifelse(source == "{PHBH_HH6%}"|target == "{PHBH_HH6%}", T, F)) %>%
                               mutate(link.with.PCL = ifelse(source == "{PCL}"|target == "{PCL}", T, F)) %>%
                               mutate(link.with.PBSA = ifelse(source == "{PBSA}"|target == "{PBSA}", T, F)) %>%
                               mutate(link.with.PBS = ifelse(source == "{PBS}"|target == "{PBS}", T, F)) %>%
                               mutate(link.with.PBAT = ifelse(source == "{PBAT}"|target == "{PBAT}", T, F)) %>%
                               mutate(link.with.monomer = ifelse(str_detect(source, "6HH|3HH|3HB|SuA|BDO|Adipate")|str_detect(target, "6HH|3HH|3HB|SuA|BDO|Adipate"), T, F)) %>%
                               mutate(link.with.monomer2 = ifelse((link.with.monomer == T & link.mic.met == T)|(link.with.monomer == T & link.with.ex == T), T, F)) %>%
                               mutate(link.with.material = ifelse(link.with.PHBH10 == T, "PHBH_HH10%",
                                                                  ifelse(link.with.PHBH6 == T, "PHBH_HH6%",
                                                                         ifelse(link.with.PCL == T, "PCL",
                                                                                ifelse(link.with.PBSA == T, "PBSA",
                                                                                       ifelse(link.with.PBS == T, "PBS",
                                                                                              ifelse(link.with.PBAT == T, "PBAT", "Other"))))))) %>%
                               arrange(G.source, G.target, source, target)
                             
                             order <-
                               c("{Test.A}", "{Test.B}", "{Test.C}", "{PCL}", "{PHBH_HH10%}", "{PHBH_HH6%}", "{PBSA}", "{PBS}", "{PBAT}", "{day.0-4}", "{day.4-14}", "{day.14-}", "{Biofilm}", "{Surface}") %>%
                               append(., res %>% filter(G.source == "Microbiome") %>% .$source %>% unique() %>% sort()) %>%
                               append(., res %>% filter(G.source == "Metabolome") %>% .$source %>% unique() %>% sort()) 
                             
                             res2 <- 
                               res %>%
                               mutate(source = factor(source, levels = order)) %>%
                               mutate(target = factor(target, levels = order)) %>%
                               arrange(source, target)
                             
                             return(res2)
                           }
                           )) %>%
  mutate(tbl_graph = map(rules_table,
                         function(x){
                           x %>%
                             as_tbl_graph() %N>%
                             mutate(id = row_number(),
                                    category = ifelse(str_detect(name, "d__"), "Microbe",
                                                      ifelse(str_detect(name, "Metabolite"), "Metabolite", "Experimental.Factor"))) %>%
                             mutate(category2 = ifelse(name == "{PHBH_HH10%}", "PHBH_HH10%",
                                                       ifelse(name == "{PHBH_HH6%}", "PHBH_HH6%",
                                                              ifelse(name == "{PCL}", "PCL",
                                                                     ifelse(name == "{PBSA}", "PBSA",
                                                                            ifelse(name == "{PBS}", "PBS",
                                                                                   ifelse(name == "{PBAT}", "PBAT",
                                                                                          ifelse(category == "Experimental.Factor", "Experimental.Factor",
                                                                                                 ifelse(category == "Metabolite", "Metabolite",
                                                                                                        ifelse(category == "Microbe", "Microbe", NA_character_)))))))))) %>%
                             mutate(category2 = factor(category2,levels = c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT", "Experimental.Factor", "Metabolite", "Microbe"))) %>%
                             mutate(is.material = ifelse(category2 %in% c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT"), "material", "other")) %>%
                             mutate(name_split = str_split(name, pattern = ";")) %>%
                             mutate(split_num = map_dbl(name_split, ~length(.))) %>%
                             mutate(lowest.ctg = map2_chr(name_split, split_num, function(x,y){x[y] %>% str_sub(., 4, -2)})) %>%
                             mutate(label = ifelse(is.material == "material", name,
                                                   ifelse(category2 == "Metabolite" & str_detect(name, "3HB|3HH|6HH|Adipate|SuA|BDO"), name, NA_character_))) %>%
                             mutate(label = str_sub(label, 2, -2)) %>%
                             mutate(label = ifelse(str_detect(label, "Metabolite_"), str_sub(label, 12, -1), label)) %>%
                             mutate(label = ifelse(str_detect(label, "Microbe_"), str_sub(label, 9, -1), label))
                         })) %>%
  mutate(graph = map(tbl_graph,
                     function(x){
                       Q <-
                         x %E>%
                         as_tibble() %>%
                         summarize(q.confidence =  quantile(confidence),
                                   q.support =  quantile(support),
                                   q.lift =  quantile(lift))
                       
                       tmp <-
                         x %E>%
                         mutate(link.with.material = factor(link.with.material, levels = c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT", "Other"))) %>%
                         filter(confidence >= Q[[1,1]] & support >= Q[[1,2]] & lift >= Q[[2,3]]) 
                         
                       #filter(link.mic.met == T) 
                       
                       link.with.material.list <-
                         tmp %E>%
                         as_tibble() %>% 
                         filter(link.with.material != "Other") %>% 
                         filter(G.source != G.target) %>%
                         gather(from, to, key = "ft", value = "id") %>% 
                         .$id %>% 
                         unique() %>%
                         as_tibble() %>%
                         set_names("id") %>%
                         mutate(link.with.material = T)
                       
                       link.with.each.material.list <-
                         tmp %E>%
                         as_tibble() %>% 
                         filter(link.with.material != "Other") %>%
                         select(from, to, link.with.PHBH10, link.with.PHBH6, link.with.PCL, link.with.PBSA, link.with.PBS, link.with.PBAT) %>%
                         gather(-from, -to, key = "link", value = "tf") %>%
                         filter(tf == TRUE) %>%
                         nest(-link) %>%
                         mutate(data = map(data, 
                                           function(x){
                                             x %>%
                                               gather(from, to, key = "ft", value = "id") %>%
                                               .$id %>%
                                               unique()
                                           })) %>%
                         unnest(data) %>%
                         rename(id = "data") %>%
                         mutate(link = str_replace(link, pattern = "link.with.", replacement = "")) %>%
                         nest(-id) %>%
                         mutate(label_text_color = map_chr(data,
                                                function(x){
                                                  str_c(x$link, collapse = "/")
                                                })) %>%
                         dplyr::select(-data)
                       
                       tmp2 <-
                         tmp %N>%
                         left_join(., link.with.material.list, by = "id") %>%
                         mutate(name2 = str_sub(name, 2, -2)) %>%
                         mutate(label2 = ifelse(link.with.material == T & category == "Metabolite", str_sub(name2, 12, -1),
                                                ifelse(link.with.material == T & category == "Microbe", lowest.ctg, NA_character_))) %>%
                         mutate(label2 = ifelse(str_detect(label2, "ROI_|unculture"), NA_character_, label2)) %>%
                         mutate(size = ifelse(is.material == "material", 1, 0.75)) %>%
                         left_join(., link.with.each.material.list, by = "id") 
                       
                       set.seed(0)
                       plot <-
                         tmp2 %N>%
                         filter(link.with.material == TRUE) %N>%
                         rename(Node = "category2") %E>%
                         rename(Link = "link.with.material") %E>%
                         filter(Link !=  "Other") %E>%
                         filter(G.source != G.target) %>%
                         ggraph(layout = "stress") +
                         geom_edge_link(color = "grey80") +
                         new_scale_color() +
                         geom_edge_link(aes(color = Link), alpha = 1) +
                         scale_edge_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                            "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                            "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                            "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                            #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                            #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16],
                                                            `NA` = "grey80")) +
                         new_scale_color() +
                         geom_node_point(aes(color = Node, size = size)) +
                         scale_size_continuous(range = c(2.5, 5)) +
                         scale_color_manual(
                           values = c(
                             "PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                             "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                             "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                             "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                             #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                             #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16],
                             #"Experimental.Factor" = brewer_pal(palette = "Greys")(9)[6], 
                             "Metabolite" = "#E0A47B", 
                             "Microbe" = "#bca89f")) +
                         new_scale_color() +
                         geom_node_text(aes(label = label2), size = 4, repel = T, check_overlap = T) +
                         theme_void() +
                         #labs(color = "Nodes") +
                         guides(size = F)
                       
                     })) %>%
  mutate(graph2 = map(tbl_graph,
                     function(x){
                       Q <-
                         x %E>%
                         as_tibble() %>%
                         summarize(q.confidence =  quantile(confidence),
                                   q.support =  quantile(support),
                                   q.lift =  quantile(lift))
                       
                       tmp <-
                         x %E>%
                         mutate(link.with.material = factor(link.with.material, levels = c("PHBH_HH10%", "PHBH_HH6%", "PCL", "PBSA", "PBS", "PBAT", "Other"))) %>%
                         filter(confidence >= Q[[1,1]] & support >= Q[[1,2]] & lift >= Q[[3,3]]) 
                       
                       #filter(link.mic.met == T) 
                       
                       link.with.material.list <-
                         tmp %E>%
                         as_tibble() %>% 
                         filter(link.with.material != "Other") %>% 
                         filter(G.source != G.target) %>%
                         gather(from, to, key = "ft", value = "id") %>% 
                         .$id %>% 
                         unique() %>%
                         as_tibble() %>%
                         set_names("id") %>%
                         mutate(link.with.material = T)
                       
                       link.with.each.material.list <-
                         tmp %E>%
                         as_tibble() %>% 
                         filter(link.with.material != "Other") %>%
                         select(from, to, link.with.PHBH10, link.with.PHBH6, link.with.PCL, link.with.PBSA, link.with.PBS, link.with.PBAT) %>%
                         gather(-from, -to, key = "link", value = "tf") %>%
                         filter(tf == TRUE) %>%
                         nest(-link) %>%
                         mutate(data = map(data, 
                                           function(x){
                                             x %>%
                                               gather(from, to, key = "ft", value = "id") %>%
                                               .$id %>%
                                               unique()
                                           })) %>%
                         unnest(data) %>%
                         rename(id = "data") %>%
                         mutate(link = str_replace(link, pattern = "link.with.", replacement = "")) %>%
                         nest(-id) %>%
                         mutate(label_text_color = map_chr(data,
                                                           function(x){
                                                             str_c(x$link, collapse = "/")
                                                           })) %>%
                         dplyr::select(-data)
                       
                       tmp2 <-
                         tmp %N>%
                         left_join(., link.with.material.list, by = "id") %>%
                         mutate(name2 = str_sub(name, 2, -2)) %>%
                         mutate(label2 = ifelse(link.with.material == T & category == "Metabolite", str_sub(name2, 12, -1),
                                                ifelse(link.with.material == T & category == "Microbe", lowest.ctg, NA_character_))) %>%
                         mutate(label2 = ifelse(str_detect(label2, "ROI_|unculture"), NA_character_, label2)) %>%
                         mutate(size = ifelse(is.material == "material", 1, 0.75)) %>%
                         left_join(., link.with.each.material.list, by = "id") 
                       
                       set.seed(0)
                       plot <-
                         tmp2 %N>%
                         filter(link.with.material == TRUE) %N>%
                         rename(Node = "category2") %E>%
                         rename(Link = "link.with.material") %E>%
                         filter(Link !=  "Other") %E>%
                         filter(G.source != G.target) %>%
                         ggraph(layout = "kk") +
                         geom_edge_link(color = "grey80") +
                         new_scale_color() +
                         geom_edge_link(aes(color = Link), alpha = 1) +
                         scale_edge_color_manual(values = c("PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                                                            "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                                                            "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                                                            "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                                                            #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                                                            #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16],
                                                            `NA` = "grey80")) +
                         new_scale_color() +
                         geom_node_point(aes(color = Node, size = size)) +
                         scale_size_continuous(range = c(2.5, 5)) +
                         scale_color_manual(
                           values = c(
                             "PHBH_HH10%" = rainbow_hcl(n=20 , c= 150, l = 40)[1], 
                             "PHBH_HH6%" = rainbow_hcl(n=20 , c= 100, l = 70)[1], 
                             "PCL" = rainbow_hcl(n=20 , c= 150, l = 70)[5], 
                             "PBSA" = rainbow_hcl(n=20 , c= 150, l = 70)[7], 
                             #"PBS" = rainbow_hcl(n=20 , c= 50, l = 70)[11], 
                             #"PBAT" = rainbow_hcl(n=20 , c= 50, l = 70)[16],
                             #"Experimental.Factor" = brewer_pal(palette = "Greys")(9)[6], 
                             "Metabolite" = "#E0A47B", 
                             "Microbe" = "#bca89f")) +
                         new_scale_color() +
                         geom_node_text(aes(label = label2), size = 4, repel = F, check_overlap = F) +
                         theme_void() +
                         #labs(color = "Nodes") +
                         guides(size = F)
                       
                     })) 

MBA_plot <- 
  MBA %>% select(graph) %>% .[[1,1]] %>% .[[1]] +
  theme(plot.margin= unit(c(1, 1, 1, 1), "lines"))  + 
  xlim(-5.5,3.5)
MBA_plot2 <- 
  MBA %>% select(graph2) %>% .[[1,1]] %>% .[[1]] +
  theme(plot.margin= unit(c(1, 1, 1, 1), "lines"))  + 
  xlim(-5.5,3.5)
ggsave(plot = MBA_plot,
       file = str_c(DIR_Integrate, "/Network.png"),
       width = 12, height = 10)


######################
### Figs for thesis
######################
ggsave(plot = MBA_plot,
       file = str_c(Dir_res_final_fig, "/Figure.3.pdf"),
       width = 9, height = 6)


  
