####Libraries####
library(Rfast)
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

# load reactome data
raw_data <- list()
filenames <- list.files("BMDExpressData/REACTOME_CSV")
for(i in 1:length(chemnames)){
  raw_data[[chemnames[i]]] <- read.csv(paste0("BMDExpressData/REACTOME_CSV/",filenames[i]), header = TRUE) %>%
    cleanupcolumns_reactome()
}

# load tPod Values
tPoD_values <- read.table(file = "tpod_values.txt")



####Filter Data####

raw_data_filtered <- list()
for(i in 1:length(raw_data)){
  raw_data_filtered[[chemnames[i]]] <- reactome_filtering(x = raw_data[[i]]) %>%
    mutate(logBMD = log10(BMD.Median)) %>%
    arrange(logBMD)
}

# add cumulative "rank" (where pathways with the same BMD get "grouped")
# This is how BMD express does it
# ranking by logBMD

for(i in 1:length(raw_data_filtered)){
  if(nrow(raw_data_filtered[[i]]>0)){
    rVal<-1
    ranker<-rVal
    cVal<-1
    counter<-cVal
    for(j in 2:nrow(raw_data_filtered[[i]])){
      if(raw_data_filtered[[i]]$logBMD[j]>raw_data_filtered[[i]]$logBMD[j-1]){
        rVal<-rVal+1+cVal-1
        ranker<-c(ranker,rVal)
        cVal<-1
        counter<-c(counter,cVal)
      }else{
        cVal<-cVal+1
        ranker<-c(ranker, rVal)
        counter<-c(counter,cVal)
      }
      
    }
    raw_data_filtered[[i]]$count<-1:nrow(raw_data_filtered[[i]])
    raw_data_filtered[[i]]$rank<-ranker
    raw_data_filtered[[i]]$rank_count<-counter
  }else{
    raw_data_filtered[[i]]$count<-vector()
    raw_data_filtered[[i]]$rank<-vector()
    raw_data_filtered[[i]]$rank_count<-vector()
  }
  
}



####Plot Data####
##Export parameters
#save tPoD figures? T or F. If F, will only display then
savefigures <- FALSE
#plot all plots together?
multi_plot <- TRUE

if(multi_plot == TRUE && savefigures == TRUE){
  png(filename = paste0("Reactome Figures/multi_plot_reactome.png"), width = 1000*length(chemnames), height = 500)
}

# x and y axis values
max_y<-max(sapply(raw_data_filtered, nrow))
min_x<- -5
max_x<- 1


for(i in 1:length(raw_data_filtered)){

plots[[i]] <- ggplot(raw_data_filtered[[i]], aes(x = logBMD, y=rank)) +
    geom_line(size=2, col="orange") +  
    geom_point(size=3, stroke=2, col="orange", shape=22, fill="white") +
    labs(title = chemnames[i], x = "logBMD", y = "Cumulative REACTOME Pathways") +
    xlim(min_x,max_x) +
    ylim(0,max_y) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = log10(as.numeric(lowestdoses[i,2])), color = "red", size = 1.5) +
    geom_vline(xintercept = tPoD_values[i, "first_mode"], color = "blue", size = 1.5) +
    geom_vline(xintercept = tPoD_values[i, "X10th_percentile"], color = "green", size = 1.5) +
    geom_vline(xintercept = tPoD_values[i, "X20th_gene"], color = "purple", size = 1.5)
}



png(filename = paste0("Reactome Figures/", "muliplot_Reactome.png"), width = 1000, height = 500)
  multiplot(plotlist=plots, cols=6)
dev.off()




#### Combined historgram and pathway accumulation plots #####


#### load BMD values for histograms ####
# load 
raw_hist_data <- list()
filenames <- list.files("BMDExpressData/BMD")

for(i in 1:length(chemnames)){
  raw_hist_data[[chemnames[i]]] <- read.table(paste0("BMDExpressData/BMD/",filenames[i]), header = TRUE, sep = "\t") %>%
    cleanupcolumns()
}

# Filter
raw_hist_filtered <- list()
for(i in chemnames){
  raw_hist_filtered[[i]] <- BMDfiltering(x = raw_hist_data[[i]], 
                                                    lowdose = lowestdoses$dose[lowestdoses$chemical==i], 
                                                    highdose = 3,
  ) %>%
    mutate(logBMD = log10(BMD))
}


# Historgram data
source("mode_antimode.R")
min_dense<-0.06  
min_bw <- 0.015 
bwFun<-"SJ"

# calculate modes (to obtain bandwidth) and breaks 
dataModes <- list()
histBreaks <- list()
for(i in chemnames){
  dataModes[[i]]<-mode.antimode(raw_hist_filtered[[i]]$logBMD ,min.size=min_dense, bw="SJ", min.bw=min_bw)
  
  histBreaks[[i]]<-seq(
    from = min(raw_hist_filtered[[i]]$logBMD) - dataModes[[i]]$bw,
    to = max(raw_hist_filtered[[i]]$logBMD) + dataModes[[i]]$bw,
    by = dataModes[[i]]$bw
  )
}

# historgram plot values
dataHist<-list()
for(i in chemnames){
  dataHist[[i]]<- hist(raw_hist_filtered[[i]]$logBMD, 
       breaks=histBreaks[[i]]
  )  
  dataHist[[i]] <- data.frame(x=dataHist[[i]]$mids, y=dataHist[[i]]$count)
}


# maximum x and y-values for multi plotting (highest freq value)
min_x <- -5
max_x <- 1

y_max <- sapply(dataHist, function(x){
    max(x$y)
  }) %>%
  max()


# Density plot values
dataDens<-list()
for(i in chemnames){
  if(!is.na(dataModes[[i]]$modes[1])){
    dataDens[[i]]<-density(raw_hist_filtered[[i]]$logBMD, bw=dataModes[[i]]$bw)
    dataDens[[i]]$y<-dataDens[[i]]$y/max(dataDens[[i]]$y)*max(dataHist[[i]]$y)
    #dataDens[[i]]<-data.frame(x=dataDens[[i]]$x, y=dataDens[[i]]$y)
  }else{
    dataDens[[i]]<-NA
  }
}




ggplot(raw_data_filtered[[i]], aes(x = logBMD, y=rank))+
  #geom_bar(data=dataHist[[i]], aes)+
  #geom_line(data=dataDens[[i]], aes(x=x, y=y), inherit.aes = FALSE, size=1.5, col=hsv(0.5,1,0.8,0.4))




  geom_line(size=2, col="orange") +  
  geom_point(size=3, stroke=2, col="orange", shape=22, fill="white") +
  labs(title = i, x = "logBMD", y = "Cumulative REACTOME Pathways") +
  xlim(min_x,max_x) +
  ylim(0,max_y) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = log10(as.numeric(lowestdoses$dose[lowestdoses$chemical==i])), color = "red", size = 1) +
  geom_vline(xintercept = tPoD_values[tPoD_values$chemical==i, "first_mode"], color = "blue", size = 1) +
  geom_vline(xintercept = tPoD_values[tPoD_values$chemical==i, "X10th_percentile"], color = "green", size = 1) +
  geom_vline(xintercept = tPoD_values[tPoD_values$chemical==i, "X20th_gene"], color = "purple", size = 1)


  
  
  
  
  i<-"DBP"
  
  ggplot(raw_hist_filtered[[i]], aes(x=logBMD)) +
    geom_histogram(binwidth = dataModes[[i]]$bw ) +
    stat_density(bw=dataModes[[i]]$bw)
  
  
  
  
