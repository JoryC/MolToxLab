#### LIBRARIES AND FUCNTION SOURCE ####
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
library(purrr)
library(ggfortify)

source("RNAseqFunctions.R")


#### IMPORT METADATA ####
metadata <- read.csv("RNAseqData/metadata2.csv", header = TRUE) %>%
  arrange(chemical, dose)


#### IMPORT SEQ DATA ####
# folders containing ONLY raw data
dataFolders <- paste0("RNAseqData/RawData/",dir("RNAseqData/RawData"))

# load data into a list
loadRaw<-list()
for(i in dataFolders){
  fileNames <- list.files(paste0(i))
  for(j in fileNames){
    loadRaw[[j]] <- read.table(paste0(i,"/",j),
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         sep = "\t",
                         strip.white = TRUE) [-c(1:4),-c(2:3)]
    colnames(loadRaw[[j]])<-c("gene", j)
  }
}

# join,
allData <- reduce(loadRaw, full_join, by="gene") %>%
# transpose,
  column_to_rownames("gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
# convert to tibble,
  as_tibble() %>%
# merge with metadata,  
  left_join(x = metadata, y = ., by="sample") %>%
# convert chemical and dose to factors,
  mutate(chemical = as.factor(chemical), dose = as.factor(dose)) %>%
# and group by chemical  
  group_by(chemical)


# NEST data: a neat little function that is similar to breaking the data into lists. However the data stay in "tibble" format, and can be manipulated with dplyr commands
nestData <- allData %>%
  nest()
  

#### FILTER ####
# uses the "countFiler' function from the "RNAseqFunctions" source code
nestData <- nestData %>%
  mutate(filterData = map(data, countFilter, grouping = "dose", median_threshold = 5)) 


#### NORMALIZE ####
nestData <- nestData %>%
  mutate(normData = map(filterData, tmmNorm)) 

# export
apply(nestData, 1, FUN = function(x){
  write_delim(x$normData,
              paste0("RNAseqData/normalizedData/",x$chemical, "_normData.txt"),
              delim = "\t"
              )
})





####PCA####

temp1 <- as.data.frame(t(norm_combined_data))
chemgroupsfactor <- as.factor(newchemgroups)
grouped_norm_data <- data.frame(chemgroupsfactor,temp1)
rm(temp1)
pca_model <- prcomp(grouped_norm_data[,-1])
autoplot(pca_model, colour = 'chemgroupsfactor', data = grouped_norm_data, label = TRUE)




