#### LIBRARIES AND FUCNTION SOURCE ####
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
# library(purrr)
library(ggfortify)
options(scipen = 9)

source("RNAseqFunctions.R")
source("BMDExpressFunctions.R") # for multiplot


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
  mutate(chemical = as.factor(chemical), dose = as.factor(dose)) %>%
# and group by chemical  
  group_by(chemical)


# NEST data: a neat little function that is similar to breaking the data into lists. However the data stay in "tibble" format, and can be manipulated with dplyr commands
nestData <- allData %>%
  nest()
  
####QC####
# Number of genes with at least 1 read per sample
ngene_per_sample <- allData %>% 
  as.data.frame() %>% 
  mutate(genecount = rowSums(.[4:ncol(allData)]!= 0)) %>%
  select(sample, chemical, dose, genecount)

#Per sample ncov5
ncov5_sample <- allData %>% 
  as.data.frame() %>% 
  mutate(genecount = rowSums(.[4:ncol(allData)] >= 5)) %>%
  select(sample, chemical, dose, genecount)

ngene_per_sample_mean <- mean(ngene_per_sample[,"genecount"])

ngene_per_chem_mean <- ngene_per_sample %>% 
  select(-sample, -dose) %>%
  group_by(chemical) %>% 
  summarise_at(vars(genecount), .funs = c(mean, sd)) %>%
  rename(ngene_mean = fn1) %>%
  rename(ngene_sd = fn2)

# Number of genes with at least 1 read across all samples
ngene_total_data <- mapply(sum, allData[,-c(1:3)])
ngene_total <- length(which(ngene_total_data != 0))

# Number of genes with a cumulative read count of at least 5 across all samples
# ncov5 <- length(which(ngene_total_data >= 5))

# Number of genes (highest to lowest read count) that make up 80% of reads
nsig80<-vector()
for(i in 1:5){
  test1 <-nestData$data[[i]] %>%
    select(-sample, -dose) %>%
    colSums()
  test1<- test1[order(test1, decreasing=TRUE)]
  test_percent <- test1/sum(test1)*100
  nsig80<-c(nsig80,sum(cumsum(test_percent)<80))
  rm(test1, test_percent)
}
names(nsig80) <- as.vector(nestData$chemical)
nsig80_mean <- mean(nsig80)

# nsig80 per sample
nsig80_sample<-vector()
for(i in 1:nrow(allData)){
  test1 <-allData[i,] %>%
    as.data.frame() %>%
    select(-sample, -dose, -chemical) %>%
    as.integer()
  test1<- test1[order(test1, decreasing=TRUE)]
  test_percent <- test1/sum(test1)*100
  nsig80_sample<-c(nsig80_sample,sum(cumsum(test_percent)<80))
  rm(test1, test_percent)
}
names(nsig80_sample) <- as.vector(allData$sample)

#Number of genes that have at least 1 read across all samples per chemical
ngene_total_chem_data <- vector()
for(i in 1:5) {
  temp <- nestData$data[[i]] %>%
    select(-sample,-dose)
  temp[temp > 0] <- 1
  temp_sum <- mapply(sum, temp)
  ngene_total_chem_data <-
    c(ngene_total_chem_data, length(which(temp_sum == nrow(nestData$data[[i]]))))
  rm(temp, temp_sum)
}

#Number of genes that have at least 1 read across all samples for all chemicals
# ngene_total_allsamples_data <- allData %>%
#   as.data.frame() %>%
#   select(-sample, -chemical, -dose)
# ngene_total_allsamples_data[ngene_total_allsamples_data > 0] <-1
# ngene_total_allsamples_sum <- mapply(sum, ngene_total_allsamples_data)
# ngene_total_allsamples <- length(which(ngene_total_allsamples_sum == nrow(allData)))

#Overall QC Summary
qcSummary <- data.frame(endpoint = c("ngene_per_sample_mean", 
                                     "ngene_total",
                                     # "ngene_total_allsamples",
                                     # "ncov5", 
                                     "nsig80_mean"),
                        value = c(ngene_per_sample_mean,
                                  ngene_total,
                                  # ngene_total_allsamples,
                                  # ncov5,
                                  nsig80_mean)
                        )

#Chemical specific QC
qcSummary_chem <- ngene_per_chem_mean %>%
  select(-ngene_sd) %>%
  mutate(ngene_allsamples = ngene_total_chem_data) %>%
  mutate(nsig80 = nsig80)

#Sample specific QC
qcSummary_sample <- ncov5_sample %>%
  rename(ncov5_sample = genecount) %>%
  mutate(nsig80_sample = nsig80_sample)
  

#### FILTER ####
# uses the "countFilter' function from the "RNAseqFunctions" source code
nestData <- nestData %>%
  mutate(filterData5 = map(data, countFilter, grouping = "dose", median_threshold = 5)) %>%
  mutate(filterData3 = map(data, countFilter, grouping = "dose", median_threshold = 3)) %>%
  mutate(filterData1 = map(data, countFilter, grouping = "dose", median_threshold = 1)) %>%
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


