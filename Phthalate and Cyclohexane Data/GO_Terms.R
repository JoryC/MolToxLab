####Libraries####
#library(Rfast)
library(tidyverse)
library(purrr)
source("BMDExpressFunctions.R")

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata_nocontrol.csv")
chemnames <- unique(metadata$chemical)
lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

####Data Import####

# NOTE: raw bmdexpress export files in .txt are broken... can be imported but not correctly and give many errors...
# Must be converted to a csv file first.

raw_data <- list()
filenames <- list.files("BMDExpressData/GO_TERM_CSV")
for(i in 1:length(chemnames)){
  raw_data[[chemnames[i]]] <- read.csv(paste0("BMDExpressData/GO_TERM_CSV/",filenames[i]), header = TRUE) %>%
    cleanupcolumns_goterm()
}

tPoD_values <- read.table(file = "tpod_values.txt")

####Filter Data####

raw_data_filtered <- list()
for(i in 1:length(raw_data)){
  raw_data_filtered[[chemnames[i]]] <- goterm_filtering(x = raw_data[[i]]) %>%
    mutate(logBMD = log10(BMD.Median)) %>%
    arrange(logBMD)
}

####Export Lowest Median GO-Term####
export_median <- TRUE

GO_Term_Endpoints <- data.frame()
for(i in 1:length(chemnames)){
  GO_Term_Endpoints[i,1] <- chemnames[i]
  GO_Term_Endpoints[i,2] <- raw_data_filtered[[i]][1,"logBMD"]
}
colnames(GO_Term_Endpoints) <- c("chemical", "go_term_median")

if(export_median == TRUE){
  write.table(GO_Term_Endpoints, file = "go_term_values.txt", quote = FALSE, sep = "\t")
}


##Plot Data####
#Export parameters
#save tPoD figures? T or F. If F, will only display then
savefigures <- TRUE
#plot all plots together?
multiplot <- FALSE

if(multiplot == TRUE && savefigures == TRUE){
  png(filename = paste0("Reactome Figures/multiplot_GO_term.png"), width = 1000*length(chemnames), height = 500)
}

if(multiplot == TRUE){
  par(mfrow = c(1,length(chemnames)))
} else {
  par(mfrow = c(1,1))
}


for(i in 1:length(raw_data_filtered)){
  if(savefigures == TRUE && multiplot == FALSE){
    png(filename = paste0("Go_Term Figures/", chemnames[i], "_GO_Term.png"), width = 1000, height = 500)
  }
  print(
    ggplot(raw_data_filtered[[i]], aes(x = logBMD)) +
      stat_ecdf(geom = "step") +
      labs(title = chemnames[i], x = "logBMD", y = "GO Term Cumulative Distribution") +
      xlim(-5,1) +
      ylim(0,1) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(xintercept = log10(as.numeric(lowestdoses[i,2])), color = "red", size = 1.5) +
      geom_vline(xintercept = tPoD_values[i, "first_mode"], color = "blue", size = 1.5) +
      geom_vline(xintercept = tPoD_values[i, "X10th_percentile"], color = "green", size = 1.5) +
      geom_vline(xintercept = tPoD_values[i, "X20th_gene"], color = "purple", size = 1.5)
  )
  if(savefigures == TRUE && multiplot == FALSE){
    dev.off()
  }
}
if(savefigures == TRUE && multiplot == TRUE){
  dev.off()
}

