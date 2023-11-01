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
#Note: Manually removed samples with poor mapping & outliers from Metadata spreadsheet & RawData Folder.

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
                             stringsAsFactors = FALSE) %>%
      rename(COUNTDATA = 5) %>%
      select(GENEID,
             BIOTYPE,
             COUNTDATA) %>%
      filter(!BIOTYPE %in% c("rRNA",
                             "Mt_rRNA")) %>%
      select(GENEID,
             COUNTDATA) %>%
      replace(is.na(.),0) %>%
      distinct()
    colnames(loadRaw[[j]])<-c("gene", j)
  }
}

loadRaw_rRNA <- list()
for(i in dataFolders){
  fileNames <- list.files(paste0(i))
  for(j in fileNames){
    loadRaw_rRNA[[j]] <- read.csv(paste0(i,"/",j),
                             header = TRUE,
                             stringsAsFactors = FALSE) %>%
      rename(COUNTDATA = 5) %>%
      select(GENEID,
             BIOTYPE,
             COUNTDATA) %>%
      filter(BIOTYPE %in% c("rRNA",
                             "Mt_rRNA")) %>%
      select(GENEID,
             COUNTDATA) %>%
      replace(is.na(.),0) %>%
      distinct()
    colnames(loadRaw_rRNA[[j]])<-c("gene", j)
  }
}

rRNA_df <- loadRaw_rRNA %>% reduce(left_join, by = "gene") %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "rRNA")) %>%
  filter(gene == "rRNA")

restofdata_df <- loadRaw %>% reduce(left_join, by = "gene") %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "not rRNA")) %>%
  filter(gene == "not rRNA")

joined_df <- full_join(rRNA_df, restofdata_df) %>%
  column_to_rownames(var="gene") %>%
  t()

write.csv(joined_df,"RNAseqData/joined_df.csv")

####

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
  left_join(x = metadata, y = ., by="sample")
# convert chemical and dose to factors,
  # mutate(chemical = as.factor(chemical), dose = as.factor(dose))
# add zeros
allData[is.na(allData)] <- 0

####TEMP####

# allData2 <- allData %>%
#   arrange(sample) %>%
#   t()
# 
# write.csv(allData2,"RNAseqData/DataCombined.csv")

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
  ylim(0, 5000) +
  labs(x = "Chemical", y = "nSig80", title = "nSig80") +
  theme_classic()

if(TRUE){
  png(file = "QC.png", width = 1000, height = 2000)
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

write.csv(QC_sample_summary, "QC/QC_sample_summary.csv", row.names = F)

######################

#create data list
data_list <- list()

for(i in as.character(unique(allData$chemical))){
  data_list[[i]] <- filter(allData, chemical == i) %>%
    select(-chemical) %>%
    column_to_rownames(var = "sample") %>%
    t()
}

##########
# R-ODAF #
##########
#Filter any samples with less than 5 million reads... no samples to filter

#Minimum gene count filter...1 gene
filter_list <- lapply(data_list, function(x) {
  x %>%
    t() %>%
    as.data.frame() %>%
    countFilter(grouping = "dose", median_threshold = 1, metadata = metadata) %>%
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
library(edgeR)

# DESeq_norm <- function(data, metadata, fact){
#   dds <- DESeqDataSetFromMatrix(countData = data,
#                                 colData = metadata,
#                                 design = ~ fact)
#   dds <- estimateSizeFactors(dds)
#   normalized_counts <- counts(dds, normalized=TRUE)
#   return(normalized_counts)
# }

norm_counts <- list()
norm_CPM_counts <- list()
all_DECounts_real <- list()

minCoverage <- 5000000 #reads
MinCount <- 1
pAdjValue <- 0.01

for(i in names(filter_list)) {
  Samp4compare <- metadata %>%
    filter(chemical == i) %>%
    filter(!dose == 0) %>%
    select(dose) %>%
    unique() %>%
    as.vector()
  names(Samp4compare) <- NULL
  Cont4compare <- metadata %>%
    filter(chemical == i) %>%
    filter(dose == 0) %>%
    select(dose) %>%
    unique() %>%
    as.vector()
  names(Cont4compare) <- NULL
  DESIGN <- "dose"
  
  sampleData <- filter_list[[i]]
  DESeqDesign <- metadata_list[[i]]
  #MAke sure these are in the exact same order!!!!!
  #all(rownames(DESeqDesign) %in% colnames(sampleData)) # Are all the row and coulm names matching?
  #all(rownames(DESeqDesign) == colnames(sampleData)) #Are they in the exact same order?
  if (all(rownames(DESeqDesign) == colnames(sampleData)) == FALSE) {
    print(
      "Metadata check: Rownames of DESeqDesign DO NOT MATCH the column names of sampleData!!!"
    )
    break
  } else {
    print("Metadata check: Rownames of DESeqDesign MATCH column names of sampleData!")
  }
  
  sampleData[is.na(sampleData)] <- 0 #Replace NA with zero
  sampleData <-
    sampleData[, (colSums(sampleData) > minCoverage)] #Exclude samples with less than 5 million reads
  
  condition1 <- Cont4compare[[1]][1]
  condition2 <- Samp4compare[[1]]
  
  DE_Design <- matrix(data = NA, ncol = 2)
  DE_Design <- DESeqDesign %>%
    filter(dose %in% c(condition1, condition2))
  DE_Design$dose <-
    factor(DE_Design$dose,
           levels = c(sort(unique(DE_Design$dose))),
           ordered = FALSE)
  samples <- sampleData[, rownames(DE_Design)]
  
  ###########
  print(paste(condition2, " vs ", condition1))
  
  colnames(samples) <- NULL
  dds <-
    DESeqDataSetFromMatrix(
      countData = round(samples),
      colData = as.data.frame(DE_Design),
      design = as.formula(paste0("~", DESIGN[1]))
    )
  
  #Re-level the data
  dds$dose <- relevel(dds$dose, ref = "0")
  
  print("Wait... (dds step executing)")
  dds <- DESeq(dds, quiet = TRUE)
  
  print(paste0(
    "Filtering genes: 75% of at least 1 group need to be above ",
    MinCount,
    " CPM"
  ))
  
  #Filter low readcounts (genes not meeting the condition : at least one condition with 75% of the samples above 1 CPM)
  
  SampPerGroup <- table(DE_Design[, DESIGN])
  Counts <- counts(dds, normalized = TRUE)
  CPMdds <- cpm(counts(dds, normalized = TRUE))
  
  Filter <- matrix(data = NA,
                   ncol = 3,
                   nrow = nrow(Counts))
  rownames(Filter) <- rownames(Counts)
  colnames(Filter) <- c("Low", "quantile", "spike")
  
  # Apply the "Relevance" condition
  
  for (gene in 1:nrow(dds)) {
    CountsPass <- NULL
    for (group in 1:length(SampPerGroup)) {
      sampleCols <-
        grep(paste0("^", dimnames(SampPerGroup)[[1]][group], "$"), DE_Design[, DESIGN]) #Once again R-ODAF has a bug where the grep is not explicit enough... for example if the group is "0" i greps all dimnames that contain a zero... annoying
      Check <-
        sum(CPMdds[gene, sampleCols] >= MinCount) >= 0.75 * SampPerGroup[group]
      CountsPass <- c(CountsPass, Check)
    }
    
    if (sum(CountsPass) > 0) {
      Filter[gene, 1] <- 1
    }	else {
      Filter[gene, 1] <- 0
    }
    
  }
  
  
  compte <- Counts[Filter[, 1] == 1,]
  Filter <- Filter[rownames(Filter) %in% rownames(compte),]
  
  print(
    paste(
      "Relevance filtering removed ",
      nrow(dds) - nrow(Filter),
      " genes from the ",
      nrow(dds),
      " assessed. ",
      nrow(Filter),
      " genes remaining",
      sep = ""
    )
  )
  
  print("Obtaining the DESeq2 results")
  
  # compute the DEGs on the genes passing the Relevance condition
  
  res_list <- NULL
  FileName_list <- NULL
  DEsamples_list <- NULL
  DECounts_list <- NULL
  Filter_DEG_list <- NULL
  for (k in condition2) {
    res <-
      results(
        dds[rownames(compte),],
        alpha = pAdjValue,
        cooksCutoff = F,
        independentFiltering = F,
        contrast = c(DESIGN[1], k, condition1),
        pAdjustMethod = 'fdr'
      )
    res_list[[paste0(i, "_", k)]] <- res
    
    FileName <- paste(i, k, "vs", condition1, "FDR", pAdjValue, sep = "_")
    FileName_list[[paste0(i, "_", k)]] <- FileName
    
    DEsamples <- subset(res, res$padj < pAdjValue)
    DEsamples_list[[paste0(i, "_", k)]] <- DEsamples
    
    DECounts <-
      compte[rownames(compte) %in% rownames(DEsamples), , drop = FALSE] #R-ODAF has a bug where the rowname is dropped if only row is returned... Fixing that with a simple argument here
    DECounts_list[[paste0(i, "_", k)]] <- DECounts
    
    Filter_DEG <-
      Filter[rownames(Filter) %in% rownames(DECounts), , drop = FALSE]
    Filter_DEG_list[[paste0(i, "_", k)]] <- Filter_DEG
  }
  
  Filter_DEGs <- do.call(rbind, Filter_DEG_list)
  rows_to_keep <- unique(rownames(Filter_DEGs))
  Filter_DEGs <- Filter_DEGs[rows_to_keep,]
  
  all_DECounts <- do.call(rbind, DECounts_list)
  rows_to_keep <- unique(rownames(all_DECounts))
  all_DECounts <- all_DECounts[rows_to_keep,]
  
  all_DEsamples <- do.call(rbind, DEsamples_list)
  all_DEsamples$comparison <- rep(names(DEsamples_list), sapply(DEsamples_list, nrow))
  
  norm_data <- counts(dds[rownames(compte)], normalized = TRUE) #Normalized counts data filtered for only relevant genes!!
  
    print("Check median against third quantile")
    print("AND")
    print("Check the presence of a spike")
    for (j in 1:length(condition2)+1) {
      for (gene in 1:nrow(all_DECounts)) {
        # Check the median against third quantile
        quantilePass <- NULL
        sampleColsg1 <-
          grep(paste0("^", dimnames(SampPerGroup)[[1]][1], "$"), DE_Design[, DESIGN])
        sampleColsg2 <-
          grep(paste0("^", dimnames(SampPerGroup)[[1]][j], "$"), DE_Design[, DESIGN])
        
        Check <- median(all_DECounts[gene,sampleColsg1]) > quantile(all_DECounts[gene,sampleColsg2], 0.75)[[1]]
        quantilePass <-c(quantilePass, Check)
        Check <- median(all_DECounts[gene,sampleColsg2]) > quantile(all_DECounts[gene,sampleColsg1], 0.75)[[1]]
        quantilePass <-c(quantilePass, Check)
        
        if ( sum(quantilePass) > 0 ) {Filter_DEGs[gene,2] <- 1 }	else { Filter_DEGs[gene,2] <- 0 }
        
        # Check for spike 
        spikePass <- NULL
        for (group in 1:length(SampPerGroup)) { 
          sampleCols<-grep(paste0("^", dimnames(SampPerGroup)[[1]][group], "$"),DE_Design[,DESIGN])
          if (max(all_DECounts[gene,sampleCols]) ==0) {Check <- FALSE} else {
            Check <- (max(all_DECounts[gene,sampleCols])/sum(all_DECounts[gene,sampleCols])) >= 1.4*(SampPerGroup[group])^(-0.66)
            spikePass<-c(spikePass, Check)
          }
        }
        if ( sum(spikePass) > 1 ) {Filter_DEGs[gene,3] <- 0 }	else { Filter_DEGs[gene,3] <- 1 }
      }
    }
    DECounts_real <- all_DEsamples[rowSums(Filter_DEGs) == 3 ,]
    DECounts_no_quant <- all_DEsamples[Filter_DEGs[,2] == 0 ,]
    DECounts_spike <- all_DEsamples[Filter_DEGs[,3] == 0 ,]
    
    print(paste("A total of ",nrow(DECounts_real), " DEGs were selected, after ",nrow(DECounts_no_quant)," genes(s) removed by the quantile rule and ", nrow(DECounts_spike)," gene(s) with a spike",sep=""))
    
    write.table(norm_data,file=paste0(here::here(), "/RNAseqData/output/", i, "_R-ODAF_DESeq2_Norm_Data.txt"), sep="\t", quote=FALSE)
    write.table(DECounts_real,file=paste0(here::here(), "/RNAseqData/output/", i, "_R-ODAF_DESeq2_DEG_table.txt"), sep="\t", quote=FALSE)
}

#Filter Gene Counts of only DEGs
#1. Import Counts data and DEG tables
norm_counts_files <- list.files(paste0(here::here(), "/RNAseqData/output"), pattern = "Norm_Data.txt")

norm_counts_list <- NULL
for (l in norm_counts_files) {
  file <- read.table(paste0(here::here(), "/RNAseqData/output/", l), row.names = 1)
  norm_counts_list[[l]] <- file
}

DEG_table_files <- list.files(paste0(here::here(), "/RNAseqData/output"), pattern = "DEG_table.txt")

DEG_table_files_list <- NULL
for (l in DEG_table_files) {
  file <- read.table(paste0(here::here(), "/RNAseqData/output/", l), row.names = 1)
  DEG_table_files_list[[l]] <- file
}

#2. Filter out only the counts rows that are in the DEG table
DEG_counts_list <- NULL
for (l in 1:length(norm_counts_files)){
  Genes <- unique(rownames(DEG_table_files_list[[l]]))
  DEG_counts <- norm_counts_list[[l]][rownames(norm_counts_list[[l]]) %in% Genes,]
  DEG_counts <- data.table::setnames(DEG_counts, old = names(DEG_counts), new = paste0(str_split(names(DEG_counts), pattern = "_", simplify = TRUE)[,4], "_", str_split(names(DEG_counts), pattern = "_", simplify = TRUE)[,5]))
  DEG_counts <- rownames_to_column(DEG_counts, var = "dose")
  DEG_counts <- data.table::setnames(DEG_counts, old = names(DEG_counts), new = str_split(names(DEG_counts), pattern = "_", simplify = TRUE)[,1]) # Be careful with this function it permanently changes colnames for norm_counts_list also. So reload in norm_counts_list if encountering errors
  DEG_counts_list[[paste0(str_split(names(norm_counts_list[l]), pattern = "_", simplify = TRUE)[,1], "_DEG_Norm_Counts_Data")]] <- DEG_counts
}

#3. Write the normalized counts, containing only DEGs (As-per the R-ODAF)
for (l in 1:length(DEG_counts_list)) {
  write.table(DEG_counts_list[[l]], paste0(here::here(), "/RNAseqData/output/", names(DEG_counts_list[l]), "R-ODAF_DESeq2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

#Rename Norm Data File columns
#1. Import Normalized Data
norm_counts_files <- list.files(paste0(here::here(), "/RNAseqData/output"), pattern = "Norm_Data.txt")
norm_counts_list <- NULL
for (l in norm_counts_files) {
  file <- read.table(paste0(here::here(), "/RNAseqData/output/", l), row.names = 1)
  norm_counts_list[[l]] <- file
}

#2. Rename
renamed_norm_counts_list <- NULL
for (l in 1:length(norm_counts_files)){
  renamed_norm_counts <- norm_counts_list[[l]]
  renamed_norm_counts <- data.table::setnames(renamed_norm_counts, old = names(renamed_norm_counts), new = paste0(str_split(names(renamed_norm_counts), pattern = "_", simplify = TRUE)[,4], "_", str_split(names(renamed_norm_counts), pattern = "_", simplify = TRUE)[,5])) #Chemical dose (4) and Replicate (5)
  renamed_norm_counts <- rownames_to_column(renamed_norm_counts, var = "dose") #Function does not work with duplicated colnames
  renamed_norm_counts <- data.table::setnames(renamed_norm_counts, old = names(renamed_norm_counts), new = str_split(names(renamed_norm_counts), pattern = "_", simplify = TRUE)[,1]) #Removing replicate
  renamed_norm_counts_list[[paste0(str_split(names(norm_counts_list[l]), pattern = "_", simplify = TRUE)[,1], "_R-ODAF_DESeq2_Norm_Data")]] <- renamed_norm_counts
  write.table(renamed_norm_counts_list[[l]], paste0(here::here(), "/RNAseqData/output/", names(renamed_norm_counts_list[l]), "_renamed_4_Tylers_WF.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}


