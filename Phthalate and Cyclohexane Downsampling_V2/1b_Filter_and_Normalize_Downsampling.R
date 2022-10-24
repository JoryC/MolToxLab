#### LIBRARIES AND FUCNTION SOURCE ####

# libraries
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
library(ggfortify)
library(scales)
library(DescTools)

#source codes
source("RNAseqFunctions.R")
source("BMDExpressFunctions.R") # for multiplot

#options
options(scipen = 9) # ensures that small numbers arent converted to scientific notation. Helpful for plotting


#### IMPORT METADATA ####
metadata <- read.csv("RNAseqData/metadata_nocontrol.csv", header = TRUE) %>%
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

####

downsample <- function(df, samplesize = NA) {
  columnnames <- colnames(df)
  temp1 <- Untable(df, freq = columnnames[2]) %>%
    unlist %>%
    unname
  temp2 <- sample(temp1, replace = T, size = samplesize)
  temp3 <- as.data.frame(table(temp2)) %>%
    setNames(., columnnames)
  rm(temp1)
  rm(temp2)
  return(temp3)
}

loadRaw_ds <- readRDS("loadRaw_ds.rds")

if(exists("loadRaw_ds") == FALSE){
  start_time <- Sys.time()
  loadRaw_ds <- lapply(loadRaw, downsample, samplesize = 12000000)
  end_time <- Sys.time()
  end_time - start_time
  saveRDS(loadRaw_ds, file = "loadRaw_ds.rds")
}


####

# join,
allData <- purrr::reduce(loadRaw_ds, full_join, by="gene") %>%
# transpose,
  column_to_rownames("gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
# convert to tibble,
  as_tibble() %>%
# merge with metadata,  
  left_join(x = metadata, y = ., by="sample")
# convert chemical and dose to factors,
  # mutate(chemical = as.factor(chemical), dose = as.factor(dose))
# add zeros
allData[is.na(allData)] <- 0

# Group by chemical and nest
nestData <- allData %>%
  group_by(chemical) %>%
  nest()

####QC####
# Total read count
readcount <- allData %>%
  as.data.frame %>%
  mutate(row_sum = rowSums(select(.,4:length(allData)))) %>%
  select(sample, chemical, dose, row_sum)

readcount_per_chem_mean <- readcount %>% 
  select(-sample, -dose) %>%
  group_by(chemical) %>% 
  summarise_at(vars(row_sum), .funs = c(mean, sd)) %>%
  rename(read_mean = fn1) %>%
  rename(read_sd = fn2)

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
  mutate(nCov5_avg=map_dbl(nCov5, ~mean(.x$nCovN))) %>%
  mutate(nCov5_sd=map_dbl(nCov5, ~sd(.x$nCovN)))
 
# nSig80: The number of the most abundant genes that make up 80% of the reads, per sample
nestData <- nestData %>%
  mutate(nSig80=map(data, ~nSig80_V2(.x, metadata = metadata))) %>%
  mutate(nSig80_avg=map_dbl(nSig80, ~mean(.x$nSig80))) %>%
  mutate(nSig80_sd=map_dbl(nSig80, ~sd(.x$nSig80)))

####Plotting####
QCplots <- list()

#readcount
QCplots[["readcount"]]<- readcount %>%
  ggplot(aes(x = chemical, y = row_sum)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
              colour = "black") +
  labs(x = "Chemical", y = "Pre-Filtered Read Count", title = "Read Count") +
  scale_y_log10(breaks = c(10^3, 10^4, 10^5, 10^6, 10^7, 10^8),
                limits = c(10^3, 10^8),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="l") +
  # geom_hline(color = "red", yintercept  = 350000) +
  theme_classic()

#number of genes per sample
QCplots[["ngene_plot"]] <-
  ngene_per_sample %>%
  ggplot(aes(x = chemical, y = genecount)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
              colour = "black") +
  ylim(0,30000)+
  labs(x = "Chemical", y = "Gene Count per Sample", title = "Gene Count") +
  theme_classic()

#nCov5 plotting
QCplots[["nCov5_plot"]] <- 
  nestData %>%
  select(chemical, nCov5) %>%
  unnest(cols=c(nCov5)) %>%
  ggplot(aes(x = chemical, y = nCovN)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
              colour = "black") +
  ylim(0, 25000) +
  # geom_hline(color = "red", yintercept  = 5000) +
  labs(x = "Chemical", y = "nCov5", title = "nCov5") +
  theme_classic()

#nsig80 plotting
QCplots[["nSig80_plot"]] <- 
  nestData %>%
  select(chemical, nSig80) %>%
  unnest(cols=c(nSig80)) %>%
  ggplot(aes(x = chemical, y = nSig80)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2,height = 0, seed = 42069),
              colour = "black") +
  # geom_hline(color = "red", yintercept  = 1000) +
  ylim(0, 4000) +
  labs(x = "Chemical", y = "nSig80", title = "nSig80") +
  theme_classic()

if(TRUE){
  png(file = "QC_ds.png", width = 500, height = 1600)
  multiplot(plotlist = QCplots, layout = matrix(c(1:4), nrow=4, byrow=TRUE))
  dev.off()
  multiplot(plotlist = QCplots, layout = matrix(c(1:4), nrow=4, byrow=TRUE))
} else {
  multiplot(plotlist = QCplots, layout = matrix(c(1:4), nrow=4, byrow=TRUE))
}



#######################################################################

QC_sample_summary <- nestData %>% 
  select(chemical, nCov5, nSig80) %>% 
  unnest(names_repair = "unique", cols = c(nCov5, nSig80)) %>%
  select(-sample...5,
         -dose...6) %>%
  rename(sample = sample...2,
         dose = dose...3) %>%
  full_join(., readcount, by = "sample") %>%
  full_join(., ngene_per_sample, by = "sample") %>%
  select(chemical, sample, dose, row_sum, genecount, nCovN, nSig80)

write.csv(QC_sample_summary, "QC/QC_sample_summary_ds.csv", row.names = F)

######################

#create data list
data_list <- list()

for(i in as.character(unique(allData$chemical))){
  data_list[[i]] <- filter(allData, chemical == i) %>%
    select(-chemical) %>%
    column_to_rownames(var = "sample") %>%
    t()
}

filter_list <- lapply(data_list, function(x) {
  x %>%
    t() %>%
    as.data.frame() %>%
    countFilter(grouping = "dose", median_threshold = 3, metadata = metadata) %>%
    select(-dose) %>%
    t()
})


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

for(i in names(filter_list)){
  dds <- DESeqDataSetFromMatrix(countData = filter_list[[i]],
                                colData = metadata_list[[i]],
                                design = ~ dose)
  dds <- estimateSizeFactors(dds)
  norm_counts[[i]] <- counts(dds, normalized=TRUE) %>%
    t()
}


for(i in names(norm_counts)){
  outData <- norm_counts[[i]] %>% as.data.frame()
  doses <- filter(metadata, chemical == i) %>%
    select(dose)
  outData <- outData %>%
    mutate(dose = doses) %>%
    relocate(dose) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column()
  write_delim(outData,
              paste0("RNAseqData/DESeqnormalizedData/",i, "_ds_DESeq_normData.txt"),
              col_names = F,
              delim = "\t")
}




# #### FILTER ####
# filterData <- nestData %>%
#   select(chemical, data) %>%
#   mutate(filterData3 = map(data, ~countFilter(.x, grouping = "dose", median_threshold = 3, metadata=metadata))) %>%
#   select(-data) %>%
#   mutate(nGene_filt3 = map_dbl(filterData3, ~nGene(.x, metadata=metadata))) %>%
#   mutate(nOverlap_filt3 = map_dbl(filterData3, ~nGeneIntercept(.x, metadata = metadata))) %>%
#   mutate(nCov5=map(filterData3, ~nCovN(.x, metadata = metadata))) %>%
#   mutate(nCov5_avg=map_dbl(nCov5, ~mean(.x$nCovN))) %>%
#   mutate(nSig80=map(filterData3, ~nSig80_V2(.x, metadata = metadata))) %>%
#   mutate(nSig80_avg=map_dbl(nSig80, ~mean(.x$nSig80)))
# 
# 
# #### NORMALIZE ####
# filterData <- filterData %>%
#   mutate(normData = map(filterData3, tmmNorm)) %>% 
#   mutate(logData = map(filterData3, logCount)) %>%
#   mutate(logNormData = map(logData, tmmNorm)) 
# 
# 
# #### EXPORT NORM DATA ####
# 
# if(TRUE){   #switch to TRUE if you want to save the output files
#   apply(filterData, 1, FUN = function(x){
#     
#     outData <- x$normData %>%
#       select(-sample, -dose) %>%
#       t() %>%
#       as.data.frame() 
#     
#     colnames(outData) <- x$normData$dose
#     outData <- data.frame(gene=rownames(outData), outData, check.names = FALSE)
#     
#     write_delim(outData,
#                 paste0("RNAseqData/normalizedData/",x$chemical, "_normData.txt"),
#                 delim = "\t"
#     )
#   })
# }


########################

# # plot before normazliation
# filterData$logData[[2]] %>%
#   select(-dose) %>%
#   column_to_rownames("sample") %>%
#   #mutate_all(~log2(.+1)) %>%
#   as.data.frame() %>%
#   t() %>%
#   boxplot()
# 
# 
# par(mfrow=c(2,3))
# for(i in 1:5){
#   # plot before normazliation
#   filterData$logNormData[[i]] %>%
#     #bind_rows() %>%
#     select(-dose) %>%
#     column_to_rownames("sample") %>%
#     #mutate_all(~log2(.+1)) %>%
#     as.data.frame() %>%
#     t() %>%
#     boxplot(horizontal=TRUE)
# }


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


