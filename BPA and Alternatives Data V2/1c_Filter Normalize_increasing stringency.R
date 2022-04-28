#### LIBRARIES AND FUCNTION SOURCE ####

# libraries
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
library(ggfortify)
library(scales)
library(DESeq2)  # just gonna load it right away, and fix the code so that it will work whether DeSeq is loaded or not

#source codes
source("RNAseqFunctions.R")
source("BMDExpressFunctions.R") # for multiplot

#options
options(scipen = 9) # ensures that small numbers arent converted to scientific notation. Helpful for plotting


#### IMPORT METADATA ####
metadata <- read.csv("RNAseqData/metadata.csv", header = TRUE) %>%
  arrange(chemical, dose)

#create metadata list for DeSeq
metadata_list <- list()

for(i in as.character(unique(metadata$chemical))){
  metadata_list[[i]] <- filter(metadata, chemical == i) %>%
    select(-chemical) %>%
    column_to_rownames(var = "sample")
}



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
allData <- purrr::reduce(loadRaw, full_join, by="gene") %>%
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

# Group by chemical and nest
nestData <- allData %>%
  group_by(chemical) %>%
  nest()

# add DeSeq metadata into tibble
nestData$DeSeqColData<-metadata_list


#### FILTERS ####

# FILTER 1: Just remove zero rows
nestData <- nestData %>%
  mutate(noZeroData = map(data, ~cleanZeroRows(.x, metadata = metadata)))

# Norm 1: }
nestData<-deSeqNorm(nestData, normCol="noZeroData", metaCol = "DeSeqColData", newCol="deSeqData")


# FILTER 2: min median of 1, NO GROUPING
nestData <- nestData %>%
  dplyr::mutate(filter2 = purrr::map(data, ~countFilter(.x, grouping = "none", metadata=metadata, median_threshold = 1)))

# Norm 2: }
nestData<-deSeqNorm(nestData, normCol="filter2", metaCol = "DeSeqColData", newCol="deSeqDataF2")




## Norm Plots
plotCol<-"deSeqDataF2"

par(mfrow=c(2,3))
for(i in 1:nrow(nestData)){
  plotData <- nestData[[plotCol]][[i]] %>%
    select(-sample, -dose)
  
  log2(plotData+min(plotData[plotData > 0])/10) %>%
    t() %>%
    graphics::boxplot(horizontal=TRUE, xlim = c(0, 20))
}



#### EXPORT NORM DATA ####

if(TRUE){   #switch to TRUE if you want to save the output files
  apply(nestData, 1, FUN = function(x){
    
    outData <- x$deSeqDataF2 %>%
      dplyr::select(-sample, -dose) %>%
      t()
    
    colnames(outData) <- x$deSeqDataF2$dose
    outData <- data.frame(gene=rownames(outData), outData, check.names = FALSE)
    
    write_delim(outData,
                paste0("RNAseqData/DESeqnormalizedData/",x$chemical, "_DeSeq_normDataF2.txt"),
                delim = "\t"
    )
  })
}


#### PCA ####
# PCA of using prcomp and plotted using autoplot

filterData <- filterData %>%
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
      labs(title = paste(chemical), colour = "Dose(µg/L)", shape = "") +
      guides(shape = FALSE) +
      theme(plot.title = element_text(hjust = 0.5))
  }))

# single plot
# nestData$PCAplots[[1]]

# multiplot
multiplot(plotlist=filterData$PCAplots, cols=3)


