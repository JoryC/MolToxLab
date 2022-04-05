#### LIBRARIES AND FUCNTION SOURCE ####

# libraries
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
library(ggfortify)

#source codes
source("RNAseqFunctions.R")
source("BMDExpressFunctions.R") # for multiplot

#options
options(scipen = 9) # ensures that small numbers arent converted to scientific notation. Helpful for plotting


#### IMPORT METADATA ####
metadata <- read.csv("RNAseqData/metadata.csv", header = TRUE) %>%
  arrange(chemical, dose)

#### IMPORT SEQ DATA ####
# folders containing ONLY raw data
dataFolders <- paste0("RNAseqData/RawData/",dir("RNAseqData/RawData"))

# load data into a list
loadRaw<-list()
for(i in dataFolders){
  fileNames <- list.files(paste0(i))
  for(j in fileNames){
    loadRaw[[j]] <- read.csv(paste0(i,"/",j),
                             header = TRUE,
                             stringsAsFactors = FALSE)
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
  mutate(chemical = as.factor(chemical), dose = as.factor(dose))

#create data list
data_list <- list()

for(i in as.character(unique(allData$chemical))){
  data_list[[i]] <- filter(allData, chemical == i) %>%
    select(-chemical,
           -dose) %>%
    column_to_rownames(var = "sample") %>%
    t()
}

# #log transform data and 
# log_data_list <- lapply(data_list, FUN = function(x){log2(x+1)})

#create metadata list
metadata_list <- list()

for(i in as.character(unique(metadata$chemical))){
  metadata_list[[i]] <- filter(metadata, chemical == i) %>%
    select(-chemical) %>%
    column_to_rownames(var = "sample")
}

library(DESeq2)

# DESeq_norm <- function(data, metadata, fact){
#   dds <- DESeqDataSetFromMatrix(countData = data,
#                                 colData = metadata,
#                                 design = ~ fact)
#   dds <- estimateSizeFactors(dds)
#   normalized_counts <- counts(dds, normalized=TRUE)
#   return(normalized_counts)
# }

norm_counts <- list()


for(i in names(data_list)){
  dds <- DESeqDataSetFromMatrix(countData = data_list[[i]],
                                colData = metadata_list[[i]],
                                design = ~ dose)
  dds <- estimateSizeFactors(dds)
  norm_counts[[i]] <- counts(dds, normalized=TRUE) %>%
    t()
}




# 
# # Group by chemical and nest
# nestData <- allData %>%
#   select(-dose) %>%
#   column_to_rownames(var = "sample") %>%
#   group_by(chemical) %>%
#   nest()
# 
# colData <- metadata %>%
#   group_by(chemical) %>%
#   nest()

# library(DESeq2)


# # Generate colData list
# colData <- list()
# for(i in unlist(unique(metadata["chemical"]))){
#   colData[[i]] <- metadata %>%
#     filter(chemical == i)
# }


