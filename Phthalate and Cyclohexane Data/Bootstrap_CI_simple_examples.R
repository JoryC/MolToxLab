#### Libraries and Source code####
library(tidyverse)
source("mode_antimode.R")

#### General Bootstrap ####
# general formal for a bootstrap is:
#   1) sample 
#   2) compute statistic
#   3) repeat


# a hypotheticl BMD dataset made by merging two "normal" data sets
logBmdValues <- c(rnorm(200, -2, 0.5), rnorm(300, 1, 0.5))
hist(logBmdValues, breaks=100)



#### Bootstrap for 20th Gene ######
set.seed(1) # (this is to ensure that the "random" sampling can be exactly replicated)
boot_20<-vector() # vector to save results
for(i in 1:2000){
  # 1:sample
    sampleData<-sample(logBmdValues, length(logBmdValues), replace=TRUE)
  # 2: compute stat (20th gene)
    s<-sampleData[order(sampleData)][20]
    boot_20 <- c(boot_20, s)
} # 3: repeat

# lower CI, Median, upper CI
quantile(boot_20, probs=c(0.025,0.5, 0.975))




#### Bootstrap for 1st Mode ######
set.seed(1) 
boot_mode<-vector()
for(i in 1:2000){
  # 1:sample
  sampleData<-sample(logBmdValues, length(logBmdValues), replace=TRUE)
  
  # 2: compute stat (1st mode)
  dataMode<-mode.antimode(sampleData ,min.size=0.06, bw="SJ", min.bw=0.015)
  s<-dataMode$modes[[1]]
  boot_mode <- c(boot_mode, s)
}

# lower CI, Median, upper CI
quantile(boot_mode, probs=c(0.025,0.5, 0.975))





#### Bootstrap for a GO Term ######

# a hypothetical list of BMDs in a GO term that has 10 genes enriched
logGObmd <- runif(n=10, min=-3, max=3)

#bootstrap
set.seed(1) 
boot_go<-vector()
for(i in 1:2000){
  # 1:sample
  sampleData<-sample(logGObmd, length(logGObmd), replace=TRUE)
  
  # 2: compute stat (median BMD of GO term)
  s<-median(sampleData)
  boot_go <- c(boot_go, s)
}

# lower CI, Median, upper CI
quantile(boot_go, probs=c(0.025,0.5, 0.975))









