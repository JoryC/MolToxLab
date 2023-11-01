####Libraries####
library(Rfast)
library(tidyverse)
library(purrr)
library(stringr)
# library(tune) #for square plots
source("BMDExpressFunctions.R")
options(scipen = 9)

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata.csv")
chemnames <- unique(metadata$chemical)
chemnames_ori <- unique(str_extract(chemnames, "[^_]+"))

lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

highestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(max)


#### REACTOME PLOT ####
# import data
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

joined_reactome <- list()
for(i in chemnames_ori){
  joined_reactome[[i]] <- inner_join(x = filtered_reactome_data[[i]], y = filtered_reactome_data[[paste0(i,"_ds")]],
                            by = "GO.Pathway.Gene.Set.Gene.ID")
}
####

joined_reactome_plot <- list()

for(i in chemnames_ori){
  if(nrow(joined_reactome[[i]]) == 0) {
    joined_reactome_plot[[i]] <- NULL
  } else {
    maxaxis <- max(c(joined_reactome[[i]]$BMD.Median.x, joined_reactome[[i]]$BMD.Median.y))
    
    joined_reactome_plot[[i]] <- joined_reactome[[i]] %>%
      ggplot(aes(x = BMD.Median.x, y = BMD.Median.y)) +
      geom_point(size = 2) +
      geom_abline(slope = 1, intercept = 0) +
      stat_summary(fun.data = mean_cl_normal) +
      geom_smooth(method = 'lm', fullrange = T) +
      xlim(0, maxaxis) +
      ylim(0, maxaxis) +
      coord_fixed(ratio = 1) +
      labs(x = "Original Median BMD", y = "Downsampled Median BMD", title = paste0(i,"(n=", nrow(joined_reactome[[i]]), " Reactome Pathways)")) +
      theme_light()
  }
}

multiplot(plotlist = joined_reactome_plot, layout = matrix(c(1:length(joined_reactome_plot)), nrow=1, byrow=TRUE))


#### GO PLOT ####
#Import Data
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
    arrange(logBMD) %>%
    as_tibble()
}

joined_go <- list()
for(i in chemnames_ori){
  joined_go[[i]] <- inner_join(x = filtered_go_data[[i]], y = filtered_go_data[[paste0(i,"_ds")]],
                                     by = "GO.Pathway.Gene.Set.Gene.ID")
}
####

joined_go_plot <- list()

for(i in chemnames_ori){
  if(nrow(joined_go[[i]]) == 0) {
    joined_go_plot[[i]] <- NULL
  } else {
    maxaxis <- max(c(joined_go[[i]]$BMD.Median.x, joined_go[[i]]$BMD.Median.y))
    
    joined_go_plot[[i]] <- joined_go[[i]] %>%
      ggplot(aes(x = BMD.Median.x, y = BMD.Median.y)) +
      geom_point(size = 2) +
      geom_abline(slope = 1, intercept = 0) +
      stat_summary(fun.data = mean_cl_normal) +
      geom_smooth(method = 'lm', fullrange = T) +
      xlim(0, maxaxis) +
      ylim(0, maxaxis) +
      coord_fixed(ratio = 1) +
      labs(x = "Original Median BMD", y = "Downsampled Median BMD", title = paste0(i,"(n=", nrow(joined_go[[i]]), " GO terms)")) +
      theme_light()
  }
}

multiplot(plotlist = joined_go_plot, layout = matrix(c(1:length(joined_go_plot)), nrow=1, byrow=TRUE))

