####Libraries####
library(tidyverse)

####nth Gene Bootstrap####
nth_gene_bootstrap <- function(x,
                               seed = 1,
                               nth_gene = 20,
                               repeats = 2000) {
  set.seed(seed)
  boot_nth_gene <- vector()
  for(i in 1:repeats){
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
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
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    s <- quantile(sampleData[order(sampleData)], probs = (nth_percent/100))
    boot_nth_percent <- c(boot_nth_percent, s)
  }
  return(quantile(boot_nth_percent, probs=c(0.025,0.5, 0.975)))
}

####First Mode Bootstrap####
mode_bootstrap <- function(x,
                           seed = 1,
                           repeats = 2000) {
  source("mode_antimode.R")
  set.seed(seed)
  boot_mode <- vector()
  for(i in 1:repeats){
    sampleData<-sample(unlist(x), length(unlist(x)), replace=TRUE)
    dataMode<-mode.antimode(sampleData, min.size=0.06, bw="SJ", min.bw=0.15)
    s<-dataMode$modes[[1]]
    boot_mode <- c(boot_mode, s)
  }
  return(quantile(boot_mode, probs=c(0.025,0.5, 0.975)))
}

####Pathway Bootstrap####
pathway_bootstrap <- function(x,
                              seed = 1,
                              repeats = 2000) {
  set.seed(seed)
  boot_pathway <- vector()
  for (i in 1:repeats) {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    s <- median(sampleData)
    boot_pathway <- c(boot_pathway, s)
  }
  return(quantile(boot_pathway, probs = c(0.025, 0.5, 0.975)))
}

####Average CI values####
averageCI <- function(x){
  lowerCIvalues <- vector()
  medianvalues <- vector()
  upperCIvalues <- vector()
  for(i in 1:length(x)){
    if(length(x) > 0){
      templow <- x[[i]][[1]]
      lowerCIvalues <- c(lowerCIvalues, templow)
      tempmed <- x[[i]][[2]]
      medianvalues <- c(medianvalues, tempmed)
      tempup <- x[[i]][[3]]
      upperCIvalues <- c(upperCIvalues, tempup)
    } else {
      lowerCIvalues = NA
      medianvalues = NA
      upperCIvalues = NA
    }
  }
  meanlowerCI <- mean(lowerCIvalues)
  meanmedianvalue <- mean(medianvalues)
  meanupperCI <- mean(upperCIvalues)
  meanvalues <- c(meanlowerCI, meanmedianvalue, meanupperCI)
  names(meanvalues) <- c("2.5%", "50%", "97.5%")
  return(meanvalues)
}



