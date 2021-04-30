setwd("E:\\Alamar Blue\\Alamar Blue R Outputs")
library(dplyr)
library(data.table)
library(abind)

#load files
baseline_1 <- read.table("E:\\Alamar Blue\\TGSH 2021_02_02\\TGSH_Baseline_1.txt", skip = 13, nrows = 8, fill = TRUE)
baseline_2 <- read.table("E:\\Alamar Blue\\TGSH 2021_02_02\\TGSH_Baseline_2.txt", skip = 13, nrows = 8, fill = TRUE)
after24h_1 <- read.table("E:\\Alamar Blue\\TGSH 2021_02_02\\TGSH_24h_1.txt", skip = 13, nrows = 8, fill = TRUE)
after24h_2 <- read.table("E:\\Alamar Blue\\TGSH 2021_02_02\\TGSH_24h_2.txt", skip = 13, nrows = 8, fill = TRUE)

#Remove last 3 columns and 7th row function
cleanup <- function(x) {
  x <- x[, -c(10:12)]
  x <- x[-7,]
  return(x)
}

#turn into a list to lapply cleanup function on all objects
templist <- list(baseline_1, baseline_2, after24h_1, after24h_2)
templist2 <- lapply(templist, cleanup)

#extract list into dataframes
baseline_1 <- data.frame(templist2[1])
baseline_2 <- data.frame(templist2[2])
after24h_1 <- data.frame(templist2[3])
after24h_2 <- data.frame(templist2[4])

#convert to array and get the mean of each cell
#along is along the "3rd dimension", i.e. not along rows or columns
baseline_avg <- rowMeans(abind(baseline_1, baseline_2, along = 3), dims=2)
after24h_avg <- rowMeans(abind(after24h_1, after24h_2, along = 3), dims=2)

#convert back to df
baseline_avg <- as.data.frame(baseline_avg)
after24h_avg <- as.data.frame(after24h_avg)

