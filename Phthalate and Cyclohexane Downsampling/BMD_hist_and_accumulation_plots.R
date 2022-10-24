####Libraries####
library(Rfast)
library(tidyverse)
library(purrr)
source("BMDExpressFunctions.R")
options(scipen = 9)

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata_nocontrol_2.csv")
chemnames <- unique(metadata$chemical)

lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

highestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(max)

#### load tPod Values ####
tPoD_values <- read.table(file = "BMDExpressData/Output/tpod_values.txt", header = T)


#### BMD HISTOGRAM #####

# load BMD values for histograms 
raw_hist_data <- raw_importer(file_path = "BMDExpressData/BMD/", 
                              file_type = ".txt", 
                              start_phrase = "Probe Id") %>%
  lapply(cleanupcolumns) %>%
  setNames(chemnames)

# Filter
filtered_hist_data <- list()
for(i in chemnames){
  filtered_hist_data[[i]] <- BMDfiltering(x = raw_hist_data[[i]], 
                                         lowdose = lowestdoses$dose[lowestdoses$chemical==i], 
                                         highdose = highestdoses$dose[highestdoses$chemical==i],
  ) %>%
    mutate(logBMD = log10(BMD))
}

# Set arguments for mode function
source("mode_antimode.R")
min_dense<-0.06  
min_bw <- 0.015 
bwFun<-"SJ"

# calculate modes (to obtain bandwidth) and breaks 
dataModes <- list()
histBreaks <- list()
for(i in chemnames){
  dataModes[[i]]<-mode.antimode(filtered_hist_data[[i]]$logBMD ,min.size=min_dense, bw="SJ", min.bw=min_bw)
  
  histBreaks[[i]]<-seq(
    from = min(filtered_hist_data[[i]]$logBMD) - dataModes[[i]]$bw,
    to = max(filtered_hist_data[[i]]$logBMD) + dataModes[[i]]$bw,
    by = dataModes[[i]]$bw
  )
}

# historgram plot values
dataHist<-list()
for(i in chemnames){
  dataHist[[i]]<- hist(filtered_hist_data[[i]]$logBMD, 
                       breaks=histBreaks[[i]],
                       plot=FALSE
  )  
  dataHist[[i]] <- data.frame(x=dataHist[[i]]$mids, y=dataHist[[i]]$count)
}

# maximum x and y-values for multi plotting (highest freq value)

# min_x <- log10(min(lowestdoses$dose))
# max_x <- log10(max(highestdoses$dose))

max_y <- sapply(dataHist, function(x){
  max(x$y)
}) %>%
  max()


# Density plot values
dataDens<-list()
for(i in chemnames){
  if(!is.na(dataModes[[i]]$modes[1])){
    dataDens[[i]]<-density(filtered_hist_data[[i]]$logBMD, bw=dataModes[[i]]$bw)
    dataDens[[i]]$y<-dataDens[[i]]$y/max(dataDens[[i]]$y)*max(dataHist[[i]]$y)
    dataDens[[i]]<-data.frame(x=dataDens[[i]]$x, y=dataDens[[i]]$y)
  }else{
    dataDens[[i]]<-NA
  }
}



histPlots<-list()
for(i in chemnames){
  histPlots[[i]]<-ggplot(filtered_hist_data[[i]], aes(x = logBMD)) +
    # x and y scales
    scale_x_continuous(limit = c(log10((lowestdoses$dose[lowestdoses$chemical == i])/10), 
                                 log10(highestdoses$dose[highestdoses$chemical == i])), oob = function(x, limits) x) + # i know this is a weird line... apparently there is a minor bug with just using "xlim" for histograms
    scale_y_continuous(limit = c(0, max_y), expand=c(0,0)) + 
    # plot data
    geom_histogram(binwidth = dataModes[[i]]$bw, color="black", fill=NA) +
    geom_line(data=dataDens[[i]], aes(x=x, y=y), inherit.aes = FALSE, size=1.5, col=hsv(0.5,1,0.8,0.4)) +
    # lines
    geom_vline(xintercept = tPoD_values[tPoD_values$chemical==i, "first_mode"], color = "red", size = 0.8, linetype = "dashed") +
    geom_vline(xintercept = tPoD_values[tPoD_values$chemical==i, "X10th_percentile"], color = "green", size = 0.8, linetype = "dotted") +
    geom_vline(xintercept = tPoD_values[tPoD_values$chemical==i, "X20th_gene"], color = "purple", size = 0.8, linetype = "dotdash") +
    # labels
    labs(title = paste0(i," (n=",nrow(filtered_hist_data[[i]]), " genes)"), x = "logBMD", y = "Count") +
    # theme
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=12))
}

multiplot(plotlist = histPlots, layout = matrix(c(1:length(histPlots)), nrow=2, byrow=TRUE))



#### REACTOME ACCUMULATION PLOT ####
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

# add cumulative "rank" (where pathways with the same BMD get the same rank)
# This is how BMD express does it
# ranking by logBMD
# NOTE after doing all this, i learned that the function "min_rank" does the same thing... ah well

for(i in 1:length(filtered_reactome_data)){
  if(nrow(filtered_reactome_data[[i]])>1){
    rVal<-1
    ranker<-rVal
    cVal<-1
    counter<-cVal
    for(j in 2:nrow(filtered_reactome_data[[i]])){
      if(filtered_reactome_data[[i]]$logBMD[j]>filtered_reactome_data[[i]]$logBMD[j-1]){
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
    filtered_reactome_data[[i]]$count<-1:nrow(filtered_reactome_data[[i]])
    filtered_reactome_data[[i]]$rank<-ranker
    filtered_reactome_data[[i]]$rank_count<-counter
  }else{
    if(nrow(filtered_reactome_data[[i]])==1){
      filtered_reactome_data[[i]]$count<-1
      filtered_reactome_data[[i]]$rank<-1
      filtered_reactome_data[[i]]$rank_count<-1
    }else{
      filtered_reactome_data[[i]]$count<-vector()
      filtered_reactome_data[[i]]$rank<-vector()
      filtered_reactome_data[[i]]$rank_count<-vector()  
    }
    
  }
  print(i)
}


# x and y axis values

# min_x <- log10(min(lowestdoses$dose))
# max_x <- log10(max(highestdoses$dose))

max_y<-max(sapply(filtered_reactome_data, nrow))



reactomePlots<-list()
for(i in chemnames){
  reactomePlots[[i]] <- ggplot(filtered_reactome_data[[i]], aes(x = logBMD, y=rank)) + # can plot by "rank" or "count"
    # x and y axis lims
    xlim(log10(lowestdoses$dose[lowestdoses$chemical == i]/10), 
         log10(highestdoses$dose[highestdoses$chemical == i])) +
    ylim(0,max_y) +
    # plot data
    geom_line(size=2, col="orange") +  
    geom_point(size=3, stroke=2, col="orange", shape=22, fill="white") +
    # labels
    labs(title = paste0(i, " (n=", nrow(filtered_reactome_data[[i]]), " pathways)"), x = "logBMD", y = "Cumulative REACTOME Pathways") +
    #theme
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=12))
}

multiplot(plotlist=reactomePlots, layout = matrix(c(1:length(reactomePlots)), nrow=2, byrow=TRUE))

# #### GO TERM ACCUMULATION PLOT ..... GROUPED BY LEVEL ####
# 
# # load GO data
# raw_go_data <- list()
# filenames <- list.files("BMDExpressData/GO_TERM_CSV")
# for(i in 1:length(chemnames)){
#   raw_go_data[[chemnames[i]]] <- read.csv(paste0("BMDExpressData/GO_TERM_CSV/",filenames[i]), header = TRUE) %>%
#     cleanupcolumns_goterm()
# }
# 
# # filter data
# filtered_go_data <- list()
# for(i in chemnames){
#   filtered_go_data[[i]] <- goterm_filtering(x = raw_go_data[[i]]) %>%
#     mutate(logBMD = log10(BMD.Median)) %>%
#     arrange(logBMD) %>%
#     as_tibble()
# }
# 
# 
# # group by GO level then ranked by logBMD
# for(i in chemnames){
#   filtered_go_data[[i]] <- filtered_go_data[[i]] %>%
#     group_by(GO.Level) %>%
#     mutate(rank=min_rank(logBMD)) %>%
#     arrange(GO.Level, logBMD) %>%
#     ungroup() %>%
#     mutate(GO.Level=as.factor(GO.Level))
# }
# 
# 
# # x and y axis values
# max_y<-max(sapply(filtered_go_data, function(x){
#   max(x$rank)
# }))
# 
# 
# goLevelPlots<-list()
# for(i in chemnames){
#   goLevelPlots[[i]] <- ggplot(filtered_go_data[[i]], aes(x = logBMD, y=rank, color=GO.Level)) + 
#     # x and y axis lims
#     xlim(log10(lowestdoses$dose[lowestdoses$chemical == i]/10), 
#          log10(highestdoses$dose[highestdoses$chemical == i])) +
#     ylim(0,max_y) +
#     # plot data
#     geom_line(size=2) +  
#     geom_point(size=3, stroke=2, shape=22, fill="white") +
#     # labels
#     labs(title = paste0(i, " (n=", nrow(filtered_go_data[[i]]), " GO terms)"), x = "logBMD", y = "Cumulative GO terms") +
#     #theme
#     theme_classic() +
#     theme(plot.title = element_text(hjust = 0.5, size=12))
# }
# 
# multiplot(plotlist=goLevelPlots, cols=3)
# 
# goLevelPlots[[2]]




#### GO TERM ACCUMULATION PLOT ..... NOT GROUPED  ####
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


# ranked by logBMD
for(i in chemnames){
  filtered_go_data[[i]] <- filtered_go_data[[i]] %>%
    mutate(rank=min_rank(logBMD)) %>%
    mutate(GO.Level=as.factor(GO.Level))
}


# x and y axis values
max_y<-max(sapply(filtered_go_data, function(x){
  max(x$rank)
}))
min_x<- -5
max_x<- 1



goPlots<-list()
for(i in chemnames){
  goPlots[[i]] <- ggplot(filtered_go_data[[i]], aes(x = logBMD, y=rank)) + 
    # x and y axis lims
    xlim(log10(lowestdoses$dose[lowestdoses$chemical == i]/10), 
         log10(highestdoses$dose[highestdoses$chemical == i])) +
    ylim(0,max_y) +
    # plot data
    geom_line(size=2, col="orange") +  
    geom_point(size=3, stroke=2, shape=22, col="orange", fill="white") +
    # labels
    labs(title = paste0(i, " (n=", nrow(filtered_go_data[[i]]), " GO terms)"), x = "logBMD", y = "Cumulative GO terms") +
    #theme
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=12))
}

multiplot(plotlist=goPlots, layout = matrix(c(1:length(goPlots)), nrow=2, byrow=TRUE))


  
