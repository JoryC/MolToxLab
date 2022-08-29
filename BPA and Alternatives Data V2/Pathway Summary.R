####Libraries####
library(tidyverse)
source("BMDExpressFunctions.R")
options(scipen = 9)


####Metadata Import####
metadata <- read.csv("RNAseqData/metadata.csv")
chemnames <- unique(metadata$chemical)

lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

highestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(max)

#### GO Import ####

raw_go_data <- raw_importer(file_path = "BMDExpressData/GO_TERM/", 
                            file_type = ".txt", 
                            start_phrase = "GO/Pathway/Gene Set/Gene ID") %>%
  lapply(cleanupcolumns_goterm) %>%
  setNames(chemnames)

# filter data
filtered_go_data <- list()
for(i in chemnames){
  filtered_go_data[[i]] <- goterm_filtering(x = raw_go_data[[i]]) %>%
    mutate(logBMD = log10(BMD.Median)) %>%
    arrange(logBMD)
}


# ranked by logBMD and filter out the top 10
for(i in chemnames){
  filtered_go_data[[i]] <- filtered_go_data[[i]] %>%
    mutate(rank=min_rank(logBMD)) %>% 
    filter(rank <= 10) %>%
    mutate(Chemical = i) %>%
    select(Chemical, everything())
}

top_10_go <- do.call(rbind, filtered_go_data) %>%
  select(Chemical,
         GO.Pathway.Gene.Set.Gene.ID,
         GO.Level,
         GO.Pathway.Gene.Set.Gene.Name,
         Fisher.s.Exact.Two.Tail,
         Percentage,
         BMD.Median,
         rank)

#### REACTOME ####

raw_reactome_data <- raw_importer(file_path = "BMDExpressData/REACTOME/", 
                            file_type = ".txt", 
                            start_phrase = "GO/Pathway/Gene Set/Gene ID") %>%
  lapply(cleanupcolumns_reactome) %>%
  setNames(chemnames)

# filter data
filtered_reactome_data <- list()
for(i in chemnames){
  filtered_reactome_data[[i]] <- reactome_filtering(x = raw_reactome_data[[i]]) %>%
    mutate(logBMD = log10(BMD.Median)) %>%
    arrange(logBMD)
}


# ranked by logBMD and filter out the top 10
for(i in chemnames){
  filtered_reactome_data[[i]] <- filtered_reactome_data[[i]] %>%
    mutate(rank=min_rank(logBMD)) %>% 
    filter(rank <= 10) %>%
    mutate(Chemical = i) %>%
    select(Chemical, everything())
}

top_10_reactome <- do.call(rbind, filtered_reactome_data) %>%
  select(Chemical,
         GO.Pathway.Gene.Set.Gene.ID,
         GO.Pathway.Gene.Set.Gene.Name,
         Fisher.s.Exact.Two.Tail,
         Percentage,
         BMD.Median,
         rank)

write.csv(top_10_go, file = "BMDExpressData/Output/go_summary.csv", row.names = F)
write.csv(top_10_reactome, file = "BMDExpressData/Output/reactome_summary.csv", row.names = F)
