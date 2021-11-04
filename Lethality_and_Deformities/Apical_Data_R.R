library(dplyr)
library(ggplot2)
### FOLDERS, FILES AND FUNCTIONS
setwd("~/MolToxLab/Lethality_and_Deformities/") #Please set your own working directory
temp1 <- list.files("Apical_Data_Sheets/", "*.csv") #List all of the files in the working directory
apicaldata <- lapply(paste0("Apical_Data_Sheets/", temp1), read.csv) #Import the data to a single object
rm(temp1)
chemnames <- gsub("\\_SR.csv$", "", list.files("Apical_Data_Sheets/", "*.csv"))
names(apicaldata) <- chemnames

####ANOVA####
anova_output <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Survival.Rate, data = apicaldata[[i]])
  anova_output[[i]] <- summary(temp2)
  rm(temp2)
}
names(anova_output) <- chemnames
print(anova_output)

####Means and SD####
apicaldatasummary <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- apicaldata[[i]] %>% group_by(Dose..mg.L.) %>%
    summarise("Survival Rate" = mean(Survival.Rate), SD = sd(Survival.Rate))
  colnames(temp2) <- c("Dose", "Survival_Rate", "SD")
  apicaldatasummary[[i]] <- as.data.frame(temp2)
  rm(temp2)
}
names(apicaldatasummary) <- chemnames
print(apicaldatasummary)


##################


for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
  ggplot(aes(x=Dose, y=Survival_Rate)) +
    geom_bar(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Survival_Rate-SD, ymax=Survival_Rate+SD), width=.2,
                  position=position_dodge(.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    ylab("Survival Rate")
    print(p)
    ggsave(names(apicaldatasummary[i]), device = "pdf", path = paste0(getwd(), "/Images"))
}
