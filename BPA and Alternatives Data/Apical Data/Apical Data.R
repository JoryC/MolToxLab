library(tidyverse)
library(DescTools)
source("Functions/Apical Functions.R")
options(scipen = 9) #prevent scientific notation... annoying while graphing
### FOLDERS, FILES AND FUNCTIONS
raw_files <- list.files(path = "RawData/", pattern = ".csv")
chemnames <- substr(raw_files, start = 1, stop = nchar(raw_files)-4)
apicaldata <- list()
for(i in chemnames){
  apicaldata[[i]] <- read.csv(paste0("RawData/", i, ".csv")) %>%
    rename(Dose = Dose..µg.L.)
}

####ANOVA####
for(i in chemnames){
  apicaldata[[i]]$Dose <- as.factor(apicaldata[[i]]$Dose)
}

anova_output <- sapply(apicaldata, combinedanova)
print(anova_output)

dunnett_output <- sapply(apicaldata, combineddunnett)
print(dunnett_output)

summarystats_output <- lapply(apicaldata, summarystats)


####Plotting####

apicalplots <- list()

for(i in chemnames){
  apicalplots[[i]] <-
    summarystats_output[[i]] %>%
    ggplot(aes(x = as.character(Dose), y = Survival_Rate, group = Dose)) +
    geom_bar(stat="identity", color ="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Survival_Rate-StdDev, ymax=Survival_Rate+StdDev), width=.2,
                  position=position_dodge(.9)) +
    theme_classic() +
    xlab("Dose (µg/L)") +
    ylab("Survival Rate") +
    scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(title = i)
}

# print(apicalplots[[1]])

multiplot(plotlist = apicalplots, cols = 2)
