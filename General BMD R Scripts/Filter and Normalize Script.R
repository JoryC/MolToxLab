#### LIBRARIES AND FUCNTION SOURCE ####
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
# library(purrr)
library(ggfortify)

source("RNAseqFunctions.R")
source("BMDExpressFunctions.R") # for multiplot


#### IMPORT METADATA ####

metadata <- read.csv("RNAseqData/metadata_nocontrol.csv", header = TRUE) %>%
  arrange(chemical, dose)

# metadata <- read.csv("RNAseqData/Control_Archive/metadata2.csv", header = TRUE) %>%
#   arrange(chemical, dose)

#### IMPORT SEQ DATA ####
# folders containing ONLY raw data

dataFolders <- paste0("RNAseqData/RawData/",dir("RNAseqData/RawData"))

# dataFolders <- paste0("RNAseqData/Control_Archive/RawData/",dir("RNAseqData/Control_Archive/RawData"))

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


#### EXPORT NORM DATA ####

if(FALSE){   #switch to TRUE if you want to save the output files
  apply(nestData, 1, FUN = function(x){
    
    outData <- x$normData %>%
      select(-sample, -dose) %>%
      t() %>%
      as.data.frame() 
    
    colnames(outData) <- x$normData$dose
    outData <- data.frame(gene=rownames(outData), outData, check.names = FALSE)
   
    write_delim(outData,
                paste0("RNAseqData/normalizedData/",x$chemical, "_normData.txt"),
                delim = "\t"
      )
  })
}


#### PCA ####
# PCA of using prcomp and plotted using autoplot

nestData <- nestData %>%
  mutate(PCAplots = map(normData, function(x){
    x <- x %>% mutate(controlcol = ifelse(dose == 0, "Control", "Treated"))
    x %>%
      select(-sample, -dose, -controlcol) %>%
      prcomp(center=TRUE, scale.=TRUE) %>%
      autoplot(data=x, 
               colour="dose", 
               # shape = "controlcol",
               size = 3) +
      theme_bw() +
      scale_shape_manual(values = c(15, 16)) +
      labs(title = paste(chemical), colour = "Dose(??M)", shape = "") +
      guides(shape = FALSE) +
      theme(plot.title = element_text(hjust = 0.5))
  }))

# single plot
# nestData$PCAplots[[1]]

# multiplot
multiplot(plotlist=nestData$PCAplots, cols=3)


