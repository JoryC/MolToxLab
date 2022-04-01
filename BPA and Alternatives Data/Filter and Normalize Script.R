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


# Group by chemical and nest
nestData <- allData %>%
  group_by(chemical) %>%
  nest()

####QC####
# Number of genes with at least 1 read for each sample
ngene_per_sample <- allData %>% 
  as.data.frame() %>% 
  mutate(genecount = rowSums(.[4:ncol(allData)]!= 0)) %>%
  select(sample, chemical, dose, genecount)

ngene_per_sample_mean <- mean(ngene_per_sample[,"genecount"])

ngene_per_chem_mean <- ngene_per_sample %>% 
  select(-sample, -dose) %>%
  group_by(chemical) %>% 
  summarise_at(vars(genecount), .funs = c(mean, sd)) %>%
  rename(ngene_mean = fn1) %>%
  rename(ngene_sd = fn2)

# Number of genes with at least 1 read across all samples
# for entire data set
allData %>%
  nGene(metadata=metadata)

# per chemical
nestData <- nestData %>%
  mutate(nGene = map_dbl(data, ~nGene(.x, metadata=metadata))) %>%
  mutate(nOverlap = map_dbl(data, ~nGeneIntercept(.x, metadata = metadata)))

# nCov5: number of genes in a sample with a count of at least 5
nestData <- nestData %>%
  mutate(nCov5=map(data, ~nCovN(.x, metadata = metadata))) %>%
  mutate(nCov5_avg=map_dbl(nCov5, ~mean(.x$nCovN)))
 
# nSig80: The number of the most abundant genes that make up 80% of the reads, per sample
nestData <- nestData %>%
  mutate(nSig80=map(data, ~nSig80_V2(.x, metadata = metadata))) %>%
  mutate(nSig80_avg=map_dbl(nSig80, ~mean(.x$nSig80)))

#### FILTER ####
filterData <- nestData %>%
  select(chemical, data) %>%
  mutate(filterData3 = map(data, ~countFilter(.x, grouping = "dose", median_threshold = 3, metadata=metadata))) %>%
  select(-data) %>%
  mutate(nGene_filt3 = map_dbl(filterData3, ~nGene(.x, metadata=metadata))) %>%
  mutate(nOverlap_filt3 = map_dbl(filterData3, ~nGeneIntercept(.x, metadata = metadata))) %>%
  mutate(nCov5=map(filterData3, ~nCovN(.x, metadata = metadata))) %>%
  mutate(nCov5_avg=map_dbl(nCov5, ~mean(.x$nCovN))) %>%
  mutate(nSig80=map(filterData3, ~nSig80_V2(.x, metadata = metadata))) %>%
  mutate(nSig80_avg=map_dbl(nSig80, ~mean(.x$nSig80)))

####Plotting####
#number of genes per sample
ngene_plot <-
  ngene_per_sample %>%
  ggplot(aes(x = chemical, y = genecount)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
              colour = "black") +
  labs(x = "Chemical", y = "Gene Count per Sample") +
  theme_classic()

#nCov5 plotting
nCov5_plot <-
  filterData$nCov5 %>%
  bind_rows() %>%
  mutate(chemical = str_extract(sample, "[^_]+")) %>%
  ggplot(aes(x = chemical, y = nCovN)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
              colour = "black") +
  ylim(0, NA) +
  labs(x = "Chemical", y = "nCov5") +
  theme_classic()

#nsig80 plotting
nSig80_plot <-
  filterData$nSig80 %>%
  bind_rows() %>%
  mutate(chemical = str_extract(sample, "[^_]+")) %>%
  ggplot(aes(x = chemical, y = nSig80)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2,height = 0, seed = 42069),
              colour = "black") +
  ylim(0, NA) +
  labs(x = "Chemical", y = "nSig80") +
  theme_classic()

#### NORMALIZE ####
filterData <- filterData %>%
  mutate(normData = map(filterData3, tmmNorm)) 



test <- filterData$normData[[2]] %>% 
  select(-names(metadata)[names(metadata) %in% names(filterData$normData[[2]])]) %>%
  unlist() %>%
  as.data.frame %>%
  rename(value = ".")

ggplot(test, aes(x = value)) +
  geom_histogram(color="black", fill=NA, bins = 1000) +
  scale_x_continuous(expand = expand_scale(0, 0))

ggplot(filterData$normData[[1]]) +
  geom_histogram(color = "black", fill = "white")


#### EXPORT NORM DATA ####

if(TRUE){   #switch to TRUE if you want to save the output files
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
      labs(title = paste(chemical), colour = "Dose(µg/L)", shape = "") +
      guides(shape = FALSE) +
      theme(plot.title = element_text(hjust = 0.5))
  }))

# single plot
# nestData$PCAplots[[1]]

# multiplot
multiplot(plotlist=nestData$PCAplots, cols=3)


