
#### FUNCTION: countFilter #####
countFilter <- function(data, grouping, median_threshold){
  require(tidyverse)
 
  medianthrs <- median_threshold

  # group and nest by "grouping"
  nestData<-data %>%
    dplyr::select(-sample) %>%
    group_by(across(all_of(grouping))) %>%
    nest()
  
  # determine median count of each nested group
  nestData <- nestData %>%
    mutate(medCount = map(data, sapply_Median))

  medianCount <- as_tibble(do.call(rbind, nestData$medCount))
  
  # minimum median count per gene
  minCount <- medianCount %>%
    sapply(min)
  
  # gene names that are above median_threshold
  keepGenes <- names(minCount)[which(minCount >= medianthrs)]
  
  # filter to only keep genes above median_threshold
  filterData <- data %>%
    dplyr::select(all_of(c("sample", "dose", keepGenes)))

  return(filterData)

  
## OLD VERSION THAT WAS SUPER SLOW BECAUSE OF "summarize_all:
##   Keeping here for reference
##   system.time test using a 19 X 32522 tibble with 6 groups
##   takes > 200 seconds!!  vs 6 seconds for "nested" method
  
# system.time({
#  medianthrs <- median_threshold
#  
#  medianCount <- data %>%
#      select(-sample) %>%
#      group_by(across(all_of(grouping))) %>%
#      summarise_all(median)
#  
#  minCount <- medianCount %>%
#    ungroup() %>%
#    select(-dose) %>%
#    summarise_all(min)
#  
#  keepGenes <- names(minCount)[which(minCount >= medianthrs)]
#  
#  filterData <- data %>%
#    select(all_of(c("sample", "dose", keepGenes)))
#  })
#  
#  
  
}


#### FUNCTION: sapply_median ####
# a function needed to combine sapply with median for nested data
sapply_Median<-function(data){
  sapply(data, median)
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
