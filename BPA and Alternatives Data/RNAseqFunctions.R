
#### FUNCTION: countFilter #####
countFilter <- function(data, grouping, median_threshold, metadata){
  require(tidyverse)
 
  # column selectors: used to select all meta data columns, so that the gene cols can be quickly selected
  metaCols <- names(metadata)[names(metadata) %in% names(data)]
  omitCols <- metaCols[!metaCols %in% grouping]
  
  # group and nest by "grouping"
  nestData<-data %>%
    dplyr::select(-all_of(omitCols)) %>%
    group_by(across(all_of(grouping))) %>%
    nest()
  
  # determine median count of each nested group
  nestData <- nestData %>%
    mutate(medCount = map(data, ~sapply(.x, median)))

  medianCount <- as_tibble(do.call(rbind, nestData$medCount))
  
  # minimum median count per gene
  minCount <- medianCount %>%
    sapply(min)
  
  # gene names that are above median_threshold
  keepGenes <- names(minCount)[which(minCount >= median_threshold)]
  
  # filter to only keep genes above median_threshold
  filterData <- data %>%
    dplyr::select(all_of(c(metaCols, keepGenes)))

  return(filterData)
}


#### Function: nCovN ####
# the number of genes in a samples with at least N counts
nCovN <- function(x, N=5,  metadata){
  metaCols <- names(metadata)[names(metadata) %in% names(x)]
  
  x %>% 
    mutate(nCovN = rowSums(.[,-which(names(x) %in% metaCols)] >= 5)) %>%
    select(all_of(metaCols), nCovN)
}

#
nCovAvg <-function(x){
  
}



#### Function: nGene ####
# counts the number of genes with a count of N in AT LEAST one sample
nGene <- function(x, metadata){  
  metaCols <- names(metadata)[names(metadata) %in% names(x)]
  
  n<- x %>%
    select(-all_of(metaCols)) %>%
    sapply(any) %>%
    sum()
  
  return(n)
}  


#### Function: nGeneIntercept ####
# counts the number of genes with a count of at least one in ALL samples
nGeneIntercept <- function(x, metadata){  
  metaCols <- names(metadata)[names(metadata) %in% names(x)]
  
  n<- x %>%
    select(-all_of(metaCols)) %>%
    sapply(all) %>%
    sum()
  
  return(n)
}  



#### FUNCTION: nSig80 ####

nSig80 <- function(x, metadata){
  metaCols <- names(metadata)[names(metadata) %in% names(x)]
  
  nsig80<-vector()
  x1 <- x %>%
    select(-all_of(metaCols))
  
  for(i in 1:nrow(x1)){
    test1 <-x1[i,]
    test1<- test1[order(test1, decreasing=TRUE)]
    test_percent <- test1/sum(test1)*100
    nsig80<-c(nsig80,sum(cumsum(test_percent)<80))
  }
  names(nsig80) <- as.vector(x$sample)
  
}

nSig80_V2 <- function(x, metadata){
  metaCols <- names(metadata)[names(metadata) %in% names(x)]
  
  x1 <- x %>%
    select(-all_of(metaCols))
  
  nSig80 <- apply(x1, 1, function(y){
    y1<- y[order(y, decreasing=TRUE)]
    y1_percent <- y1/sum(y1)*100
    n<-sum(cumsum(y1_percent)<80)
    return(n)
  })
  
  output<-x %>%
    select(all_of(metaCols)) %>%
    cbind(nSig80)
}




#### FUNCTION: tmmNorm ####
tmmNorm <-function(data){
  require(edgeR)
  require(tidyverse)
  
  # transpose and convert to tibble
  testData <- data %>%
    dplyr::select(-dose,) %>%
    column_to_rownames("sample") %>%
    t() %>%
    as.data.frame()

  # normalize
  sizeFactor <- calcNormFactors(testData, method = "TMM")
  normData<-testData
  for(i in 1:ncol(normData)){
    normData[,i]<-normData[,i]/(sum(normData[,i])*sizeFactor[i])
  }  

  # transpose back to orignal format
  normData <- normData %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    as_tibble()
  
  # map back to sample and dose
  newData <- data %>%
    dplyr::select(sample, dose) %>%
    left_join(normData, by="sample") 
  
  return(newData)
  
}
