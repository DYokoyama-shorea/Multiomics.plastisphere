rm(list = ls())

### parameter
min.read <- 5000   # minimum total read per sample for filtering samples with low quality  


### package
library(tidyverse)
library(readxl)
library(vegan)

###

CWD <- getwd()
dir_qiime_out <- str_c(CWD, "/data/qiime_out", sep = "")
dir_data_out <-  str_c(CWD, "/data", sep = "")

DF_orig <-
  as_tibble_col(dir_qiime_out, column_name = "dir") %>%
  mutate(info = map(dir,
                    function(x){
                      read_xlsx(str_c(x, "/metadata_plastisphere.xlsx"), sheet = "all") %>%
                        mutate(Description = as.character(Description))
                    })) %>%
  mutate(data_asv = map(dir, 
                        function(x){
                          read_tsv(str_c(x, "/asv_table.tsv", sep = ""), skip = 1)
                        })) %>%
  mutate(taxonomy = map(dir, 
                        function(x){
                          read_tsv(str_c(x, "/taxonomy.tsv", sep = ""))
                        })) %>%
  mutate(data_pathway = map(dir, 
                            function(x){
                              read_tsv(str_c(x, "/pathway_table.tsv", sep = ""), skip = 1)
                            }))  %>%
  mutate(data_ko = map(dir, 
                       function(x){
                         read_tsv(str_c(x, "/ko_table.tsv", sep = ""), skip = 1)
                       }))  %>%
  mutate(data_ec = map(dir,
                       function(x){
                         read_tsv(str_c(x, "/ec_table.tsv", sep = ""), skip = 1)
                       }))  #%>%
  #mutate(data_prokatlas = map(dir,
  #                            function(x){
  #                              file <- dir_ls(str_c(x, "/prokatlas"), regexp = "prokatlas.txt") 
  #                              read_tsv(file)
  #                            }))  


DF <-
  DF_orig %>%
  mutate(sample.for.analysis = map(data_asv,
                                   function(x){
                                     colSums(x[-1]) %>%
                                       as_tibble(rownames = "sample") %>%
                                       filter(value > min.read) %>%
                                       .$sample
                                   })) %>%
  mutate(data_asv.2 = pmap(list(data_asv, info, sample.for.analysis), 
                          function(x, y, z){
                            data.x <- x %>% rename(ID = `#OTU ID`) %>% select(one_of(c("ID", z)))
                            info <-
                              y %>%
                              rename(sample = `#SampleID`) %>%
                              select(sample, Project, Description)
                            res <-
                              data.x %>%
                              rename(Taxon = "ID") %>%
                              select(Taxon, everything()) %>%
                              gather(-Taxon, key = "sample", value = "read") %>%
                              group_by(sample) %>%
                              mutate(sum.read = sum(read)) %>%
                              ungroup() %>%
                              mutate(prop = read / sum.read) %>%
                              left_join(., info, by ="sample")
                          })) %>%
  mutate(data_asv_rarefied = map2(data_asv, sample.for.analysis,
                                 function(x, y){
                                   tmp <-
                                     x %>%
                                     rename(ID = `#OTU ID`) %>% 
                                     select(one_of(c("ID", y))) %>%
                                     gather(-ID, key = "sample", value = "read") %>%
                                     spread(key = ID, value = read)
                                   rarefy.limit <- 
                                     tmp %>%
                                     select(-sample) %>%
                                     rowSums() %>% 
                                     as_tibble() %>% 
                                     filter(value > min.read) %>% 
                                     arrange(value) %>% 
                                     .[[1,1]]
                                     
                                   res <-
                                     tibble(seed = seq(1,10,1)) %>%
                                     mutate(data.rarefied = map(seed,
                                                                function(xx){
                                                                  set.seed(xx)
                                                                  
                                                                  rrarefy(tmp %>% select(-sample), rarefy.limit) %>%
                                                                    as_tibble() %>%
                                                                    bind_cols(tmp %>% select(sample), .)
                                                                }))
                                   return(res)
                                 })) %>%
  mutate(data_asv_taxon = pmap(list(data_asv, taxonomy, info, sample.for.analysis),
                               function(x, y, z, a){
                                 data.x <- x %>% rename(ID = `#OTU ID`) %>% select(one_of(c("ID", a)))
                                 data.y <- y %>% rename(ID = `Feature ID`)
                                 info <-
                                   z %>%
                                   rename(sample = `#SampleID`) %>%
                                   select(sample, Project, Description)
                                 
                                 res <-
                                   data.x %>%
                                   left_join(., data.y, by = "ID") %>%
                                   select(-ID, -Confidence) %>%
                                   select(Taxon, everything()) %>%
                                   gather(-Taxon, key = "sample", value = "read") %>%
                                   filter(!str_detect(Taxon, "d__Eukaryota"))  %>%             # delete Eukayota
                                   group_by(Taxon, sample) %>%
                                   summarize(read = sum(read)) %>%
                                   ungroup() %>%
                                   group_by(sample) %>%
                                   mutate(sum.read = sum(read)) %>%
                                   ungroup() %>%
                                   mutate(prop = read / sum.read) %>%
                                   left_join(., info, by ="sample")
                                   
                               })) %>%
  mutate(data_pathway.2 = pmap(list(data_pathway, sample.for.analysis, info),
                               function(x, y, z){
                                 info <-
                                   z %>%
                                   rename(sample = `#SampleID`) %>%
                                   select(sample, Project, Description)
                                 tmp1 <-
                                   x %>%
                                   rename(ID = `#OTU ID`) %>%
                                   select(one_of("ID", y)) %>%
                                   gather(-ID, key = "sample", value = "read")
                                 tmp2 <-
                                   tmp1 %>%
                                   group_by(sample) %>%
                                   summarize(sum.read = sum(read))
                                 res <-
                                   left_join(tmp1, tmp2, by = "sample") %>%
                                   mutate(prop = read/ sum.read) %>%
                                   left_join(., info, by ="sample")
                                 
                                 return(res) 
                               })) %>%
  mutate(data_ko.2 = pmap(list(data_ko, sample.for.analysis, info),
                          function(x, y, z){
                            info <-
                              z %>%
                              rename(sample = `#SampleID`) %>%
                              select(sample, Project, Description)
                            tmp1 <-
                              x %>%
                              rename(ID = `#OTU ID`) %>%
                              select(one_of("ID", y)) %>%
                              gather(-ID, key = "sample", value = "read")
                            tmp2 <-
                              tmp1 %>%
                              group_by(sample) %>%
                              summarize(sum.read = sum(read))
                            res <-
                              left_join(tmp1, tmp2, by = "sample") %>%
                              mutate(prop = read/ sum.read) %>%
                              left_join(., info, by ="sample")
                            
                            return(res) 
                          })) %>%
  mutate(data_ec.2 = pmap(list(data_ec, sample.for.analysis, info),
                          function(x, y, z){
                            info <-
                              z %>%
                              rename(sample = `#SampleID`) %>%
                              select(sample, Project, Description)
                            tmp1 <-
                              x %>%
                              rename(ID = `#OTU ID`) %>%
                              select(one_of("ID", y)) %>%
                              gather(-ID, key = "sample", value = "read")
                            tmp2 <-
                              tmp1 %>%
                              group_by(sample) %>%
                              summarize(sum.read = sum(read))
                            res <-
                              left_join(tmp1, tmp2, by = "sample") %>%
                              mutate(prop = read/ sum.read) %>%
                              left_join(., info, by ="sample")
                            
                            return(res) 
                          })) #%>%
  #mutate(data_prokatlas.2 = pmap(list(data_prokatlas, sample.for.analysis, info),
  #                               function(x, y, z){
  #                                 info <-
  #                                   z %>%
  #                                   rename(sample = `#SampleID`) %>%
  #                                   select(sample, Project, Description)
  #                                 tmp1 <-
  #                                   x %>%
  #                                   rename(ID = "X1") %>%
  #                                   select(one_of("ID", y)) %>%
  #                                   gather(-ID, key = "sample", value = "read")
  #                                 tmp2 <-
  #                                   tmp1 %>%
  #                                   group_by(sample) %>%
  #                                   summarize(sum.read = sum(read))
  #                                 res <-
  #                                   left_join(tmp1, tmp2, by = "sample") %>%
  #                                   mutate(prop = read/ sum.read) %>%
  #                                   left_join(., info, by ="sample")
  #                               }))


#################################################################################################################
###
###
### taxonomy
###
###
#################################################################################################################

asv.table_tmp <-
  DF %>%
  select(dir, data_asv_taxon) %>%
  unnest() %>%
  select(dir, Taxon, sample, Project, Description, read, sum.read) 

DF2 <-  
  tibble(category = c("d", "p", "c", "o", "f", "g", "s"),
         num = seq(1,7,1),
         data = rep(nest(asv.table_tmp)[[1]],7)) %>%
  mutate(data2 = map2(data, num,
                      function(x, y){
                        x %>%
                          mutate(taxa = str_split(Taxon, pattern = ";")) %>%
                          mutate(taxa = map(taxa, function(x){x[1:y]})) %>%
                          #                          mutate(taxa = map(taxa, function(x){x[!is.na(x)]})) %>% # to summarize unassigned higher levels
                          mutate(taxa = map_chr(taxa, function(x){str_c(x, collapse = ";")})) %>%
                          mutate(taxa = ifelse(is.na(taxa), "Unassigned", taxa)) %>%
                          group_by(taxa, dir, sample, Project, Description, sum.read) %>%
                          summarize(read = sum(read)) %>%
                          mutate(prop = read/sum.read) %>%
                          ungroup()
                      })) %>%
  mutate(data3 = map(data2,
                     function(x){
                       x %>%
                         select(-sum.read, -read) %>%
                         spread(key = taxa, value = prop) %>%
                         mutate_all(~ifelse(is.na(.), 0, .))
                     })) %>%
  mutate(data4 = map(data2,
                     function(x){
                       x %>%
                         select(-sum.read, -read) %>%
                         spread(key = taxa, value = prop) %>%
                         mutate_all(~ifelse(is.na(.), 0, .))
                     })) %>%
  select(category, data3)




DF2_asv <-  
  DF %>% 
  select(dir, data_asv.2) %>%
  unnest() %>%
  select(dir, Taxon, sample, Project, Description, read, sum.read)  %>%
  nest() %>%
  mutate() %>%
  mutate(category = "asv", num = 8) %>%
  mutate(data2 = map2(data, num,
                      function(x, y){
                        x %>%
                          mutate(taxa = Taxon) %>%
                          select(-Taxon) %>%
                          #group_by(taxa, dir, sample, Project, Description, sum.read) #%>%
                          #summarize(read = sum(read)) %>%
                          mutate(prop = read/sum.read) %>%
                          ungroup() %>%
                          select(taxa, dir, sample, Project, Description, sum.read, read, prop)
                      })) %>%
  mutate(data3 = map(data2,
                     function(x){
                       x %>%
                         select(-sum.read, -read) %>%
                         spread(key = taxa, value = prop) #%>%
                         #mutate_all(~ifelse(is.na(.), 0, .))
                     })) %>%
  select(category, data3)


DF3 <-
  bind_rows(DF2, DF2_asv)

DF_rarefied <- 
  DF %>% 
  select(info, taxonomy, sample.for.analysis, data_asv_rarefied) %>%
  unnest(data_asv_rarefied) %>%
  mutate(data_asv_taxon_rarefied = pmap(list(data.rarefied, taxonomy, info),
                                        function(x, y, z){
                                          data.y <- y %>% rename(ID = `Feature ID`)
                                          info <-
                                            z %>%
                                            rename(sample = `#SampleID`) %>%
                                            select(sample, Project, Description)
                                          
                                          res <-
                                            x %>%
                                            gather(-sample, key = "ID", value = read) %>%
                                            left_join(., data.y, by = "ID") %>%
                                            select(-ID, -Confidence) %>%
                                            select(Taxon, everything()) %>%
                                            filter(!str_detect(Taxon, "d__Eukaryota"))  %>%             # delete Eukayota
                                            group_by(Taxon, sample) %>%
                                            summarize(read = sum(read)) %>%
                                            ungroup() %>%
                                            group_by(sample) %>%
                                            mutate(sum.read = sum(read)) %>%
                                            ungroup() %>%
                                            mutate(prop = read / sum.read) %>%
                                            left_join(., info, by ="sample")
                                        })) 
  
  



write_rds(DF3, str_c(dir_data_out, "/microbe.rds"))
write_rds(DF_rarefied, str_c(dir_data_out, "/microbe_rarefied.rds"))

#################################################################################################################
###
###
### PICRUSt
###
###
#################################################################################################################

DF_PICRUSt <-
  DF %>%
  select(dir, data_pathway.2, data_ko.2, data_ec.2) %>%
  gather(data_pathway.2, data_ko.2, data_ec.2, key = "db", value = "data") %>%
  unnest() %>%
  select(-read, -sum.read) %>%
  nest(-db) %>%
  mutate(data2 = map(data, 
                     function(x){
                       x %>%
                         spread(key = ID, value = prop)
                     }))

write_rds(DF_PICRUSt %>% select(-data), str_c(dir_data_out, "/DF_PICRUSt.rds"))
