####Libraries####
library(tidyverse)
source("BMDExpressFunctions.R")

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata.csv")
chemnames <- unique(metadata$chemical)

lowestdoses <- unique(metadata[, c("chemical", "dose")]) %>%
  group_by(chemical) %>%
  filter(dose > 0) %>%
  summarise_all(min)

highestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(max)

####GO Term BMD Extraction####
#Import
go_raw_data <- raw_importer(file_path = "BMDExpressData/GO_TERM/", 
                                           file_type = ".txt", 
                                           start_phrase = "GO/Pathway/Gene Set/Gene ID") %>%
  lapply(cleanupcolumns_goterm) %>%
  setNames(chemnames)


#Filtering
go_raw_data_filtered <- list()
for (i in 1:length(go_raw_data)) {
  go_raw_data_filtered[[chemnames[i]]] <-
    goterm_filtering(x = go_raw_data[[i]]) %>%
    mutate(logBMD = log10(BMD.Median)) %>%
    arrange(logBMD)
}

#Filtering 2: Electric Boogaloo
go_BMD_list <- list()
for (i in chemnames) {
  go_BMD_list[[i]] <- BMD_list_extraction(go_raw_data_filtered[[i]])
}

go_BMD_list_logtransformed <- list()
for(i in chemnames){
  go_BMD_list_logtransformed[[i]] <- lapply(go_BMD_list[[i]], log10)
}

####Reactome BMD Extraction####
#Import
reactome_raw_data <- raw_importer(file_path = "BMDExpressData/REACTOME/", 
                                                    file_type = ".txt", 
                                                    start_phrase = "GO/Pathway/Gene Set/Gene ID") %>%
  lapply(cleanupcolumns_reactome) %>%
  setNames(chemnames)

#Filtering
reactome_raw_data_filtered <- list()
for (i in 1:length(reactome_raw_data)) {
  reactome_raw_data_filtered[[chemnames[i]]] <-
    reactome_filtering(x = reactome_raw_data[[i]]) %>%
    mutate(logBMD = log10(BMD.Median)) %>%
    arrange(logBMD)
}

#Filtering 2: Electric Boogaloo
reactome_BMD_list <- list()
for (i in chemnames) {
  reactome_BMD_list[[i]] <-
    BMD_list_extraction(reactome_raw_data_filtered[[i]])
}

reactome_BMD_list_logtransformed <- list()
for(i in chemnames){
  reactome_BMD_list_logtransformed[[i]] <- lapply(reactome_BMD_list[[i]], log10)
}

####Full BMD list####

BMD_raw_data <- raw_importer(file_path = "BMDExpressData/BMD/", 
                             file_type = ".txt", 
                             start_phrase = "Probe Id") %>%
  lapply(cleanupcolumns) %>%
  setNames(chemnames)

# Filter
all_BMD_list_logtransformed <- list()
for (i in chemnames) {
  all_BMD_list_logtransformed[[i]] <-
    BMDfiltering(x = BMD_raw_data[[i]],
                 lowdose = lowestdoses$dose[lowestdoses$chemical == i],
                 highdose = highestdoses$dose[highestdoses$chemical == i]) %>%
    mutate(logBMD = log10(BMD)) %>%
    select(logBMD) %>%
    pull() %>%
    list()
}

####Export Data#### 
saveRDS(go_BMD_list_logtransformed, "BMDExpressData/RDS/go_BMD_list_logtransformed.RDS")
saveRDS(reactome_BMD_list_logtransformed, "BMDExpressData/RDS/reactome_BMD_list_logtransformed.RDS")
saveRDS(all_BMD_list_logtransformed, "BMDExpressData/RDS/all_BMD_list_logtransformed.RDS")

