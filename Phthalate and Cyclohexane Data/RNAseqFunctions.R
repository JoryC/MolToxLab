#### FUNCTION: countFilter #####
countFilter <- function(data, grouping, median_threshold){
  require(tidyverse)
  
  medianthrs <- median_threshold
  
  medianCount<-data %>%
    group_by(across(all_of(grouping))) %>%
    select(-sample) %>%
    mutate_all(as.numeric) %>%   #NOTE: if left as "integer" it is SUPER SLOW
    summarize_all(median)
  
  minCount <- medianCount %>%
    ungroup() %>%
    select(-dose) %>%
    summarise_all(min)
  
  keepGenes <- names(minCount)[which(minCount >= medianthrs)]
  
  filterData <- data %>%
    select(all_of(c("sample", "dose", keepGenes)))
  
  return(filterData)
}