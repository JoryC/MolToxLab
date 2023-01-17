#### LOAD LIBRARIES ####
library(tidyverse)
library(dplyr)
library(abind)
library(car)
library(ggplot2)
library(data.table)
library(readr)
library(tibble)

#### LOAD FILE AND CLEAN UP ####
#NOTE: run getwd() to confirm you are in the ~/MolToxLab/Alamar_Blue directory
# names of all files in all folders (one folder for now)
folderNames<-list.files("Data/")
# start on first folder for now (get it working for TGSH for now, 5th folder)
folder<-folderNames
files_in_folder <- list.files(paste0("Data/",folder))

#### simulation of loop on all chemicals ####
finalResults<-vector()
for(i in 1:length(folderNames)){
  p_results<-c(0.07, 0.09)
  names(p_results)<-c("dose","replicate")
  finalResults<-rbind(finalResults, p_results)
}
row.names(finalResults)<-folderNames


#### Apply Tyler's (altered by Jory) script to import all chemicals in 2 dfs ####

#Subset files
B_1 <- grep("*Baseline_1.txt",files_in_folder)
B_2 <- grep("*Baseline_2.txt",files_in_folder)
h24_1 <- grep("*24h_1.txt",files_in_folder)
h24_2 <- grep("*24h_2.txt",files_in_folder)

#Load Files
baseline_1 <- lapply(X = paste0("Data/",folder,"/", files_in_folder[B_1]), FUN = read.table, skip = 13, nrows = 8, fill = TRUE)
baseline_2 <- lapply(X = paste0("Data/",folder,"/", files_in_folder[B_2]), FUN = read.table, skip = 13, nrows = 8, fill = TRUE)
after24h_1 <- lapply(X = paste0("Data/",folder,"/", files_in_folder[h24_1]), FUN = read.table, skip = 13, nrows = 8, fill = TRUE)
after24h_2 <- lapply(X = paste0("Data/",folder,"/", files_in_folder[h24_2]), FUN = read.table, skip = 13, nrows = 8, fill = TRUE)

#Remove last 3 columns and 7th row function
cleanup <- function(x) {
  x <- x[, -c(10:12)]
  x <- x[-7,]
  return(x)
}

#turn into a list to lapply cleanup function on all objects
templist <- list(baseline_1, baseline_2, after24h_1, after24h_2)
baseline_1_list <- lapply(X = templist[[1]], FUN = cleanup)
baseline_2_list <- lapply(X = templist[[2]], FUN = cleanup)
h24_1_list <- lapply(X = templist[[3]], FUN = cleanup)
h24_2_list <- lapply(X = templist[[4]], FUN = cleanup)

#extract list into dataframes
baseline_1 <- rbindlist(baseline_1_list, idcol = FALSE)
baseline_2 <- rbindlist(baseline_2_list, idcol = FALSE)
after24h_1 <- rbindlist(h24_1_list, idcol = FALSE)
after24h_2 <- rbindlist(h24_2_list, idcol = FALSE)

rm(templist, baseline_1_list, baseline_2_list, h24_1_list, h24_2_list) #remove templists

#convert to array and get the mean of each cell
#along is along the "3rd dimension", i.e. not along rows or columns
baseline_avg <- rowMeans(abind(baseline_1, baseline_2, along = 3), dims=2)
after24h_avg <- rowMeans(abind(after24h_1, after24h_2, along = 3), dims=2)

#convert back to df
baseline_avg <- as.data.frame(baseline_avg)
after24h_avg <- as.data.frame(after24h_avg)

rm(after24h_1, after24h_2, baseline_1, baseline_2) #remove old separated data

#seperate out data from control (i.e. final row) NOTE:Breaks if run multiple times...
baseline_avg_control <- baseline_avg[c(7,14,21,28,35),1:4]
baseline_avg_control <- t(baseline_avg_control)
baseline_avg <- t(baseline_avg[-c(7,14,21,28,35),])
after24h_avg_control <- after24h_avg[c(7,14,21,28,35),1:4]
after24h_avg_control <- t(after24h_avg_control)
after24h_avg <- t(after24h_avg[-c(7,14,21,28,35),])

#assign highest dose and reformatting data #

#dose assignment
options(scipen = 999) #Make sure R doesn't convert to scientific notation
alldoses <- as.vector(read_tsv("Highest_Dose.txt", col_names = TRUE, skip = 1,)) #Located in Alamar_Blue directory

#Repeat Each iteration of alldoses 9 times and create a vector

Doses_column <- vector(mode = "double")
for (i in 1:nrow(alldoses)) {
  for (j in 1:ncol(alldoses)) {
    x <- rep(alldoses[i, j], each = 9)
    Doses_column <- as.vector(c(Doses_column, x), mode = "double")
  }
}

#Formatting Data into columns function - edited
formatting <- function(x) {
  value <- as.vector((x)) #transpose data into a vector
  temp1 <- data.frame(value) #convert into a dataframe
  #add doses
  temp1$dose <- Doses_column
  #add reps, rep function repeats string, "each" is how often each value is repeated before going to the next one
  temp1$replicate <- rep(c("A", "B", "C"), each = 3)
  return(temp1)
}

baseline <- formatting(baseline_avg)
after24h <- formatting(after24h_avg)

#### NORMALIZE ####

#test data
Delta <- after24h$value - baseline$value
temp1 <- as.data.frame(cbind(Delta))
temp1$dose <- Doses_column 
temp1$replicate <- rep(c("A", "B", "C"), each = 3)
temp1$group <- rep(c("BPA", "BPAF", "DES", "EE2", "TGSH"), each = 54)
rm(Delta)

#### ANCOVA ####

#set factor levels in ascending order
temp1$dose <- ordered(temp1$dose,
                         levels = rev(unique(temp1$dose)))

#summary of data
summary1 <- group_by(temp1, group) %>% #(dataframe, group)
  summarise(
    count = n(),
    mean = mean(Delta, na.rm = TRUE), #replace deltavalue with data, na.rm is NA skipping
    sd = sd(Delta, na.rm = TRUE) #replace deltavalue with data, na.rm is NA skipping
  )

#check for homogeneity of variance
output <- split(temp1, temp1$group)
LeveneResults <- t(as.data.frame(c(
  leveneTest(y = Delta ~ dose, data = output$BPA),
  leveneTest(y = Delta ~ dose, data = output$BPAF),
  leveneTest(y = Delta ~ dose, data = output$DES),
  leveneTest(y = Delta ~ dose, data = output$EE2),
  leveneTest(y = Delta ~ dose, data = output$TGSH)
)
  )
) 
Insig <- print(all(LeveneResults[c(2,5,8,11,14),1] > 0.05))
#Insig should = TRUE

#ANCOVA
#run ancova, replicate is covariable, should not use type 1 since multiple factors
bpaancova <- aov(formula = Delta ~ dose + replicate, data = output$BPA)
summary(bpaancova)
bpafancova <- aov(formula = Delta ~ dose + replicate, data = output$BPAF)
summary(bpafancova)
desancova <- aov(formula = Delta ~ dose + replicate, data = output$DES)
summary(desancova)
ee2ancova <- aov(formula = Delta ~ dose + replicate, data = output$EE2)
summary(ee2ancova)
tgshancova <- aov(formula = Delta ~ dose + replicate, data = output$TGSH)
summary(tgshancova)
#aov(values ~ groups + covariable, data = x)

#ANOVA
bpaanova <-Anova(bpaancova, type = "3")
summary(bpaanova)
bpafanova <- Anova(bpafancova, type = "3")
summary(bpafanova)
desanova <- Anova(desancova, type = "3")#intercept being significant just means grand mean =/= 0
summary(desanova)
ee2anova <- Anova(ee2ancova, type = "3")
summary(ee2anova)
tgshanova <- Anova(tgshancova, type = "3")
summary(tgshanova)

#### PLOTTING ####
# (still testing)
Final <- temp1 %>%
  mutate(mean = rep((summary1$mean), each = 54)) %>%
  mutate(sd = rep((summary1$sd), each = 54))

ggplot() +
  geom_boxplot(data = subset(Final, group == "BPA"), aes(x = dose, y = Delta)) +
  geom_boxplot(data = subset(Final, group == "BPAF"), aes(x = dose, y = Delta)) +
  geom_boxplot(data = subset(Final, group == "DES"), aes(x = dose, y = Delta)) +
  geom_boxplot(data = subset(Final, group == "EE2"), aes(x = dose, y = Delta)) +
  geom_boxplot(data = subset(Final, group == "TGSH"), aes(x = dose, y = Delta)) +
  facet_wrap( ~ group, scales = "free")
