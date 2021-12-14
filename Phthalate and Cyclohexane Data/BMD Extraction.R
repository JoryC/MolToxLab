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

raw_data_filtered <- list()
for(i in 1:length(raw_data)){
  raw_data_filtered[[chemnames[i]]] <- BMDfiltering(x = raw_data[[i]], 
                                 BMD.div.BMDL = 20, 
                                 BMDU.div.BMDL = 40, 
                                 BMDU.div.BMD = 20, 
                                 lowdose = lowestdoses$dose[i], 
                                 highdose = 3,
                                 fitP=0.1
                                 )
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


