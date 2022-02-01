####Libraries####
library(tidyverse)
source("BMDExpressFunctions.R")

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata_nocontrol.csv")
chemnames <- unique(metadata$chemical)
lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

filenames <- list.files(pattern = ".*_values.txt")
columnnames <- c("chemical", "GO Term Median", "Reactome Median", "First Mode", "10th Percentile", "20th Gene")

####Create summary of data####
BMD_summary <- data.frame()
for(i in 1:length(chemnames)){
  BMD_summary[i,1]<- chemnames[i]
}
colnames(BMD_summary) <- "chemical"

for(i in 1:length(filenames)){
  temp1 <-  read.table(filenames[i])
  BMD_summary <- merge(BMD_summary, temp1, by = 'chemical')
  rm(temp1)
}

colnames(BMD_summary) <- columnnames

BMD_summary_temp <- BMD_summary %>% 
  column_to_rownames(var = 'chemical') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column()

BMD_Summary_List <- list()
for(i in 1:length(chemnames)){
  BMD_Summary_List[[chemnames[i]]] <- BMD_summary_temp[,c(1,i+1)]
  colnames(BMD_Summary_List[[chemnames[i]]]) <- c("Endpoint_Type", "logBMD")
}


min_x <- -5
max_x <- 1

summary_plots <- list()

for(i in chemnames){
  summary_plots[[i]] <- ggplot(BMD_Summary_List[[i]], aes(x = logBMD, y = Endpoint_Type)) + # can plot by "rank" or "count"
    # x and y axis lims
    xlim(min_x,max_x) +
    # plot data
    geom_point(size=3, col="Black") +
    # labels
    labs(title = paste(names(BMD_Summary_List[i]), "Endpoints"), x = "logBMD", y = "") +
    # theme
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=12))
}

multiplot(plotlist=summary_plots, cols=3)
