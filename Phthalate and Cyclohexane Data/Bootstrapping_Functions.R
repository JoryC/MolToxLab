####Libraries####
library(tidyverse)
source("mode_antimode.R")

####nth Gene Bootstrap####
nth_gene_bootstrap <- function(x,
                               seed = 1,
                               nth_gene = 20,
                               repeats = 2000) {
  set.seed(seed)
  boot_nth_gene <- vector()
  for(i in 1:repeats){
    sampleData <- sample(x, length(x), replace = TRUE)
    s <- sampleData[order(sampleData)][nth_gene]
    boot_nth_gene <- c(boot_nth_gene, s)
  }
  return(quantile(boot_nth_gene, probs=c(0.025,0.5, 0.975)))
}


####nth Percentile Bootstrap####
nth_percent_bootstrap <- function(x,
                                  seed = 1,
                                  nth_percent = 10,
                                  repeats = 2000) {
  
  set.seed(seed)
  boot_nth_percent <- vector()
  for(i in 1:repeats){
    sampleData <- sample(x, length(x), replace = TRUE)
    s <- quantile(sampleData[order(sampleData)], probs = (nth_percent/100))
    boot_nth_percent <- c(boot_nth_percent, s)
  }
  return(quantile(boot_nth_percent, probs=c(0.025,0.5, 0.975)))
}

####First Mode Bootstrap####
mode_bootstrap <- function(x,
                           seed = 1,
                           repeats = 2000,
                           minsize = 0.6,
                           
                           ) {
  
  set.seed(seed)
  boot_mode <- vector()
}








logBMDvalues <- readRDS("all_BMD_list_logtransformed.RDS")

metadata <- read.csv("RNAseqData/metadata_nocontrol.csv")
chemnames <- unique(metadata$chemical)
lowestdoses <- unique(metadata[, c("chemical", "dose")]) %>%
  group_by(chemical) %>%
  filter(dose > 0) %>%
  summarise_all(min)

testlist <- list()

for(i in chemnames){
  testlist[[i]] <- nth_percent_bootstrap(unlist(logBMDvalues[[i]]))
}