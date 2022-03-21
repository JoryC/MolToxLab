library(dplyr)
library(ggplot2)
### FOLDERS, FILES AND FUNCTIONS
# setwd("~/MolToxLab/Lethality_and_Deformities/") # Please set your own working directory
temp1 <- list.files("Apical_Data_Sheets/", "*.csv") #List all of the files in the working directory
apicaldata <- lapply(paste0("Apical_Data_Sheets/", temp1), read.csv) #Import the data to a single object
rm(temp1)
chemnames <- gsub("\\_SR.csv$", "", list.files("Apical_Data_Sheets/", "*.csv"))
names(apicaldata) <- chemnames

####ANOVA####
#Survival Rate
anova_output_SR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Survival.Rate, data = apicaldata[[i]])
  anova_output_SR[[i]] <- summary(temp2)
  rm(temp2)
}
names(anova_output_SR) <- chemnames
print(anova_output_SR)
#Deformity Rate
anova_output_DR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Deformity.Rate, data = apicaldata[[i]])
  anova_output_DR[[i]] <- summary(temp2)
  rm(temp2)
}
names(anova_output_DR) <- chemnames
print(anova_output_DR)
#Hatch Rate
anova_output_HR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ apicaldata[[i]][[5]], data = apicaldata[[i]])
  anova_output_HR[[i]] <- summary(temp2)
  rm(temp2)
}
names(anova_output_HR) <- chemnames
print(anova_output_HR)

####Means and SD####
apicaldatasummary <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- apicaldata[[i]] %>% group_by(Dose..mg.L.) %>%
    summarise("Survival Rate" = mean(Survival.Rate), SD_SR = sd(Survival.Rate),
              "Deformity Rate" = mean(Deformity.Rate), SD_DR = sd(Deformity.Rate),
              "Hatch Rate" = mean(apicaldata[[i]][[5]]), SD_HR = sd(apicaldata[[i]][[5]]))
  colnames(temp2) <- c("Dose", "Survival_Rate", "SR_SD", "Deformity_Rate", "DR_SD", "Hatch_Rate", "HR_SD")
  apicaldatasummary[[i]] <- as.data.frame(temp2)
  rm(temp2)
}
names(apicaldatasummary) <- chemnames
print(apicaldatasummary)


##################

#Survival Rate
#Saving Images
for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
  ggplot(aes(x=Dose, y=Survival_Rate)) +
    geom_bar(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Survival_Rate-SR_SD, ymax=Survival_Rate+SR_SD), width=.2,
                  position=position_dodge(.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    ylab("Survival Rate")
    print(p)
    ggsave(paste0(names(apicaldatasummary[i]), ".png"), device = "png", path = paste0(getwd(), "/Images/Survival_Rate"))
}

#Deformity Rate
for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
    ggplot(aes(x=Dose, y=Deformity_Rate)) +
    geom_bar(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    ylab("Rate of Deformity")
  print(p)
  ggsave(paste0(names(apicaldatasummary[i]), ".png"), device = "png", path = paste0(getwd(), "/Images/Deformity_Rate"))
}

#Hatch Rate
for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
    ggplot(aes(x=Dose, y=Hatch_Rate)) +
    geom_bar(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Hatch_Rate-HR_SD, ymax=Hatch_Rate+HR_SD), width=.2,
                  position=position_dodge(.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    ylab("Hatch Rate")
  print(p)
  ggsave(names(apicaldatasummary[i]), device = "pdf", path = paste0(getwd(), "/Images/Hatch_Rate"))
}
