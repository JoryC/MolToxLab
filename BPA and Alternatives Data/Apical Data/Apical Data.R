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

lethal_anova_output <- sapply(apicaldata, combinedanova, type = "lethal")
print(lethal_anova_output)

deform_anova_output <- sapply(apicaldata, combinedanova, type = "deform")
print(deform_anova_output)

lethal_dunnett_output <- sapply(apicaldata, combineddunnett, type = "lethal")
print(lethal_dunnett_output)

deform_dunnett_output<- sapply(apicaldata, combineddunnett, type = "deform")
print(deform_dunnett_output)

summarystats_output <- lapply(apicaldata, summarystats)


####Plotting####
#lethal plots
lethal_plots <- list()

# for(i in chemnames){
#   lethal_plots[[i]] <-
#     summarystats_output[[i]] %>%
#     ggplot(aes(x = as.character(Dose), y = Survival_Rate, group = Dose)) +
#     geom_bar(stat="identity", color ="black", fill="white",
#              position=position_dodge()) +
#     geom_errorbar(aes(ymin=Survival_Rate-Survival_StdDev, ymax=Survival_Rate+Survival_StdDev), width=.2,
#                   position=position_dodge(.9)) +
#     geom_jitter(aes(x = as.character(Dose), y = apicaldata[[i]]$Survival.Rate),
#                 position = position_jitter(width = 0.2, height = 0, seed = 42069),
#                 colour = "black") +
#     theme_classic() +
#     xlab("Dose (µg/L)") +
#     ylab("Survival Rate") +
#     scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#     labs(title = i)
# }

# for(i in chemnames) {
#   lethal_plots[[i]] <-
#     apicaldata[[i]] %>%
#     ggplot(aes(x = as.character(Dose), y = Survival.Rate, group = Dose)) +
#     geom_boxplot(outlier.shape = NA, width = 0.5) +
#     geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 2022),
#                 colour = "black") +
#     theme_classic() +
#     xlab("Dose (µg/L)") +
#     ylab("Survival Rate") +
#     scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#     labs(title = i)
# }

for(i in chemnames) {
  lethal_plots[[i]] <-
    apicaldata[[i]] %>%
    ggplot(aes(x = as.character(Dose), y = Survival.Rate, group = Dose)) +
    stat_summary(fun = mean, geom = "bar", 
                 color = "black", 
                 fill = "white", 
                 position=position_dodge(), 
                 width = 0.8) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - sd(x),
                 fun.max = function(x) mean(x) + sd(x),
                 geom = "errorbar",
                 width=.2, 
                 position=position_dodge(.9)) +
    geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 2022),
                colour = "black") +
    theme_classic() +
    xlab("Dose (µg/L)") +
    ylab("Survival Rate") +
    scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(title = i)
}

multiplot(plotlist = lethal_plots, cols = 2)

#deform plots
deform_plots <- list()

# for(i in chemnames){
#   deform_plots[[i]] <-
#     summarystats_output[[i]] %>%
#     ggplot(aes(x = as.character(Dose), y = Deform_Rate, group = Dose)) +
#     geom_bar(stat="identity", color ="black", fill="white", 
#              position=position_dodge()) +
#     geom_errorbar(aes(ymin=Deform_Rate-Deform_StdDev, ymax=Deform_Rate+Deform_StdDev), width=.2,
#                   position=position_dodge(.9)) +
#     geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
#                 colour = "black") +
#     theme_classic() +
#     xlab("Dose (µg/L)") +
#     ylab("Deformity Rate") +
#     scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#     labs(title = i)
# }

# for(i in chemnames) {
#   deform_plots[[i]] <-
#     apicaldata[[i]] %>%
#     ggplot(aes(x = as.character(Dose), y = Deformity.Rate, group = Dose)) +
#     geom_boxplot(outlier.shape = NA, width = 0.5) +
#     geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 42069),
#                 colour = "black") +
#     theme_classic() +
#     xlab("Dose (µg/L)") +
#     ylab("Deformity Rate") +
#     scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#     labs(title = i)
# }


for(i in chemnames) {
  deform_plots[[i]] <-
    apicaldata[[i]] %>%
    ggplot(aes(x = as.character(Dose), y = Deformity.Rate, group = Dose)) +
    stat_summary(fun = mean, geom = "bar", 
                 color = "black", 
                 fill = "white", 
                 position=position_dodge(), 
                 width = 0.8) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - sd(x),
                 fun.max = function(x) mean(x) + sd(x),
                 geom = "errorbar",
                 width=.2, 
                 position=position_dodge(.9)) +
    geom_jitter(position = position_jitter(width = 0.2, height = 0, seed = 2022),
                colour = "black") +
    theme_classic() +
    xlab("Dose (µg/L)") +
    ylab("Deformity Rate") +
    scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(title = i)
}


multiplot(plotlist = deform_plots, cols = 2)
