####Libraries####
library(tidyverse)
library(purrr)
source("functions.R")

####Metadata Import####
metadata <- read.csv("metadata.csv")
colnames(metadata) <- c("Chemical","0.3","0.03","0.003","0.0003","0.00003","3","0")
chemnames <- metadata[1:6,1]
lowestdoses <- c(0.0003,0.0003,0.0003,0.0003,0.00003,0.0003) #need better way to do this...

####Data Import####
raw_data <- list()
filenames <- list.files(pattern = "BMD_")

for(i in 1:length(chemnames)){
  raw_data[[chemnames[i]]] <- read.table(filenames[i], fill = TRUE, header = TRUE, sep = "\t")
}

raw_data <- lapply(raw_data, cleanupcolumns)
raw_data_filtered <- list()

for(i in 1:length(raw_data)){
  raw_data_filtered[[chemnames[i]]] <- BMDfiltering(x = raw_data[[i]], 
                                 BMD.div.BMDL = 20, 
                                 BMDU.div.BMDL = 40, 
                                 BMDU.div.BMD = 20, 
                                 lowdose = lowestdoses[i], 
                                 highdose = 3)
}


# raw_data_filtered_2 <- lapply(
#   raw_data,
#   FUN = BMDfiltering,
#   BMD.div.BMDL = 20,
#   BMDU.div.BMDL = 40,
#   BMDU.div.BMD = 20,
#   lowdose = lowestdoses,
#   highdose = 3
# )
# 
# names(raw_data_filtered_3) = chemnames
# 
# raw_data_filtered_2 <- lapply(chemnames, FUN = function(y){
#   BMDfiltering(x=raw_data[[y]], 
#                BMD.div.BMDL = 20,
#                BMDU.div.BMDL = 40,
#                BMDU.div.BMD = 20,
#                lowdose = lowestdoses[y],
#                highdose = 3)
#   })


BMDs_filtered <- lapply(raw_data_filtered, select, BMD)




combined_BMDs <- data.frame()

for(i in 1:length(raw_data_filtered)){
  combined_BMDs <- data.frame(combined_BMDs, raw_data_filtered[[1]] %>% select(BMD))
}


