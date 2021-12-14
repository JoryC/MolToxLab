####Libraries####
library(tidyverse)
library(purrr)
source("BMDExpressFunctions.R")

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata2.csv")
chemnames <- unique(metadata$chemical)
lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

####Data Import####
raw_data <- list()
filenames <- list.files("BMDExpressData/BMD")

for(i in 1:length(chemnames)){
  raw_data[[chemnames[i]]] <- read.table(paste0("BMDExpressData/BMD/",filenames[i]), header = TRUE, sep = "\t") %>%
    cleanupcolumns()
}

#### FILTER DATA ####

raw_data_filtered <- list()
for(i in 1:length(raw_data)){
  raw_data_filtered[[chemnames[i]]] <- BMDfiltering(x = raw_data[[i]], 
                                 BMD.div.BMDL = 20, 
                                 BMDU.div.BMDL = 40, 
                                 BMDU.div.BMD = 20, 
                                 lowdose = lowestdoses$dose[i], 
                                 highdose = 3,
                                 fitP=0.1
                                 ) %>%
    mutate(logBMD = log10(BMD))
}



#### MODE ANALYSIS #####


# load source code
source("mode_antimode.R")


# variables and options (current values based on optimiztion results)
min_dense<-0.055  # minimum probability density to be considered a "mode"
min_bw <- 0.015   # minimum bandwidth (too much resolution gives strange "peaks")
bwFun<-"SJ" # choose nrd0 or SJ. the "bandwidth" function to use to determine modes. I've selected the Sheather & Jones (1991) method. see: https://www.ncbi.nlm.nih.gov/pubmed/24885339

# lists tosave results
dataModes <- list()
firstMode <- list()
firstAntiMode <- list()
histBreaks <- list()


# run in a loop!
for(i in 1:length(raw_data_filtered)){
  
  # calculate modes
  dataModes[[i]]<-mode.antimode(raw_data_filtered[[i]]$logBMD ,min.size=min_dense, bw="SJ", min.bw=min_bw)
  
  # First Mode
  firstMode[[i]] <- dataModes[[i]]$modes[1]
  
  # First Antimode
  firstAntiMode[[i]] <-dataModes[[i]]$anti.modes[1]
  
  # histogram plot with modes/antimodes using the same bandwidth as the mode algorithm
  histBreaks[[i]]<-seq(
    from = min(raw_data_filtered[[i]]$logBMD) - dataModes[[i]]$bw,
    to = max(raw_data_filtered[[i]]$logBMD) + dataModes[[i]]$bw,
    by = dataModes[[i]]$bw
  )
  
  
}

for(i in 1:length(raw_data_filtered)){
  #i<-1
  hist(raw_data_filtered[[i]]$logBMD, breaks=histBreaks[[i]], prob=TRUE, main=names(raw_data_filtered)[i])
  lines(density(raw_data_filtered[[i]]$logBMD, bw=dataModes[[i]]$bw), col=hsv(0.5,1,0.8,0.4), lwd=3)
  abline(v=dataModes[[i]]$modes, col="red", lwd=3)
  abline(v=dataModes[[i]]$anti.modes, col="Blue", lwd=3)
}













