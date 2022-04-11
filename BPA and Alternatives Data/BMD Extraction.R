####Libraries####
library(Rfast)
library(tidyverse)
library(purrr)
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


####Data Import####
raw_data <- list()
filenames <- list.files("BMDExpressData/BMD")

for (i in 1:length(chemnames)) {
  raw_data[[chemnames[i]]] <-
    read.table(paste0("BMDExpressData/BMD/", filenames[i]),
               header = TRUE,
               sep = "\t") %>%
    cleanupcolumns()
}

#### FILTER DATA ####

raw_data_filtered <- list()
for(i in 1:length(raw_data)){
  raw_data_filtered[[chemnames[i]]] <- BMDfiltering(x = raw_data[[i]], 
                                 lowdose = lowestdoses$dose[i], 
                                 highdose = highestdoses$dose[i],
                                 ) %>%
    mutate(logBMD = log10(BMD))
}


#### MODE ANALYSIS #####

# load source code
source("mode_antimode.R")


# variables and options (current values based on optimization results)
min_dense<-0.06  # minimum probability density to be considered a "mode"
min_bw <- 0.015 # minimum bandwidth (too much resolution gives strange "peaks")
bwFun<-"SJ" # choose nrd0 or SJ. the "bandwidth" function to use to determine modes. I've selected the Sheather & Jones (1991) method. see: https://www.ncbi.nlm.nih.gov/pubmed/24885339
# does bwFun work here or down at below at "# calculate modes"

# lists tosave results
dataModes <- list()
firstMode <- list()
firstAntiMode <- list()
histBreaks <- list()


# run in a loop!
for(i in 1:length(raw_data_filtered)){
  
  # calculate modes
  dataModes[[i]]<-mode.antimode(raw_data_filtered[[i]]$logBMD ,min.size=min_dense, bw="SJ", min.bw=min_bw)
  
  # # First Mode
  # firstMode[[i]] <- dataModes[[i]]$modes[1]
  # 
  # # First Antimode
  # firstAntiMode[[i]] <-dataModes[[i]]$anti.modes[1]
  
  # histogram plot with modes/antimodes using the same bandwidth as the mode algorithm
  histBreaks[[i]]<-seq(
    from = min(raw_data_filtered[[i]]$logBMD) - dataModes[[i]]$bw,
    to = max(raw_data_filtered[[i]]$logBMD) + dataModes[[i]]$bw,
    by = dataModes[[i]]$bw
  )
}

tpod_values <- data.frame(unique(metadata[,"chemical"])) %>% 
  rename(chemical = contains("chemical"))

for(i in 1:length(raw_data)){
  
  #first mode tPoD
  tpod_values[i,"first_mode"] <- dataModes[[i]]$modes[1]
  
  #calculate 10th percentile tPoD
  tpod_values[i,"10th_percentile"] <- raw_data_filtered[[i]][,"BMD"] %>% 
    quantile(probs = 0.1) %>% log10()
  
  #calculate 20th gene tPoD
  if(length(raw_data_filtered[[i]][,"BMD"]) >= 20){
    tpod_values[i,"20th_gene"] <- raw_data_filtered[[i]][,"BMD"] %>%
      Rfast::nth(20, descending = FALSE) %>% log10()
  } else {
    tpod_values[i,"20th_gene"] <- NA
  }
}

####Data Export####

#export tPoDs
# write.table(tpod_values, file = "BMDExpressData/Output/tpod_values.txt", quote = FALSE, sep = "\t")


##Export parameters
#save tPoD figures? T or F. If F, will only display then
savefigures <- F
#plot all plots together?
multiplot <- TRUE

if(multiplot == TRUE && savefigures == TRUE){
  png(filename = paste0("tPoD Figures/multiplot_tPoD.png"), width = 1000*length(chemnames), height = 500)
}

if(multiplot == TRUE){
  par(mfrow = c(1,length(chemnames)))
} else {
  par(mfrow = c(1,1))
}

for(i in 1:length(raw_data_filtered)){
  #i<-1
  if(savefigures == TRUE && multiplot == FALSE){
    png(filename = paste0("tPoD Figures/", chemnames[i], "_tPoD.png"), width = 1000, height = 500)
  }
  hist(raw_data_filtered[[i]]$logBMD,
       breaks=histBreaks[[i]],
       prob=TRUE,
       main=names(raw_data_filtered)[i],
       xlim = c(log10(lowestdoses$dose[i]),log10(highestdoses$dose[i])),
       xlab = "logBMD"
       )
  lines(density(raw_data_filtered[[i]]$logBMD, bw=dataModes[[i]]$bw), col=hsv(0.5,1,0.8,0.4), lwd=3)
  abline(v=log10(lowestdoses[i,"dose"]), col="red", lwd = 3)
  abline(v=tpod_values[i, "first_mode"], col = "blue", lwd = 3)
  abline(v=tpod_values[i, "10th_percentile"], col = "green", lwd = 3)
  abline(v=tpod_values[i, "20th_gene"], col = "purple", lwd = 3)
  if(savefigures == TRUE && multiplot == FALSE){
    dev.off()
  }
}
if(savefigures == TRUE && multiplot == TRUE){
  dev.off()
}


