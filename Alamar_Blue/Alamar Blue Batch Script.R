####Libraries####
library(data.table)
library(abind)
library(dplyr)
library(tidyr)
library(readr)
library(car)
library(purrr)
library(broom)
library(naniar)
library(ggplot2)

####Directory####
getwd() #Output should be "/*/*/MoltToxLab/Alamar_Blue" or something similar
folderNames <- list.files("Data/Sub/")#Character list of all folders in /Data/Sub/
folder <- folderNames
files_in_folder <- list.files(paste0("Data/Sub/", folder)) #Character list of all files in each folder (alphabetical order)

####Import the Data####
list <- list() #Create an empty list
for (k in 1:length(files_in_folder)) {
    list[[k]] <- fread(file = paste0("Data/All/", files_in_folder[[k]]), 
                            header = FALSE,
                            skip = 13,
                            nrows = 8,
                            select = c(3:11)
    )
}
names(list) <- files_in_folder #Give each dataframe/datatable its corresponding name
#Create a quick functioin to remove the pesky empty row
cleanup <- function(x) {
  x <- x[-7,]
  return(x)
}
#Remove those pesky rows
list <- lapply(list, FUN = cleanup)
#Subset data
#Data is listed in alphabetical order in 'list' object
#REVIEW: Not ideal way to subset data... naming rows by Chemical - Dose can fix this and make subsetting more streamlined
Baseline_1 <- c(seq(from = 3, to = 104, by = 4)) #WARNING: seq() to argument is dynamic... will change with more data
Baseline_2 <- c(seq(from = 4, to = 104, by = 4))
h24_1 <- c(seq(from = 1, to = 104, by = 4)) #leading h because object can't start with numeric
h24_2 <- c(seq(from = 2, to = 104, by = 4))
#Create your 4 dataframes to prepare for creating average data sets
Base_1_dfl <- rbindlist(list[Baseline_1])
Base_2_dfl <- rbindlist(list[Baseline_2])
h24_1_dfl <- rbindlist(list[h24_1])
h24_2_dfl <- rbindlist(list[h24_2])
#Simple average of the two data sets (Day 0 - Baseline, and Day 1 - 24h)
Baseline_avg <- as.data.frame(
  rowMeans(
    abind(Base_1_dfl, Base_2_dfl, along = 3), #Along '3rd dimension'... not row or column
    dims=2 
    )
  ) 
h24_avg <- as.data.frame(
  rowMeans(
    abind(h24_1_dfl, h24_2_dfl, along = 3), 
    dims=2
    )
  )
rm(Base_1_dfl, Base_2_dfl, h24_1_dfl, h24_2_dfl, Baseline_1, Baseline_2, h24_1, h24_2, list)
#Subset the control data (First 4 cells by column of every 7th row)
Base_control <- Baseline_avg[c(seq(from = 7, to = 182, by = 7)),c(1:4)] #WARNING: seq() to argument is dynamic... will change with more data
row.names(Base_control) <- folderNames
h24_control <- h24_avg[c(seq(from = 7, to = 182, by = 7)),c(1:4)]
row.names(h24_control) <- folderNames
#Remove the control row of the *_avg dfs
Baseline_avg <- Baseline_avg[-c(seq(from = 7, to = 182, by = 7)),] #WARNING
h24_avg <- h24_avg[-c(seq(from = 7, to = 182, by = 7)),]
#Calculate the Delta of the Day 0 vs Day 1
Delta <- h24_avg - Baseline_avg
colnames(Delta) <- c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3", "C_1", "C_2", "C_3")

#Assign dose names to all corresponding chemicals
Delta_1 <- Delta %>%
  mutate(Dose = rep(c("Dose_1", "Dose_2", "Dose_3", "Dose_4", "Dose_5", "Dose_6"), times = 26), 
         Chemical = rep(folderNames, times = 1, each  = 6)) %>%
  group_by(Chemical) %>%
  select(Chemical, everything())

#Tidy up the data
Delta_2 <- Delta_1 %>%
  gather(key = "Replicate", value = "Fluorescence", A_1:C_3, factor_key = TRUE)  %>%
  separate(Replicate, into = c("Group", "Replicate"), convert = TRUE) %>%
  arrange(Chemical)

#Assign Dose Values
AllDoses <- read.csv(file = "Highest_Dose_Aug.csv", skip = 1,  header = TRUE)
AllDoses_2 <- AllDoses %>%
  gather(key = Dose, value = "Dose(mg/L)", Dose_1:Dose_6)

#Final Data Frame
Tidy_Data <- Delta_2 %>%
  inner_join(AllDoses_2) %>%
  select(Chemical, Group, Replicate, 'Dose(mg/L)', Fluorescence)
#Separate Date From Chemical Name
Tidy_Data <- Tidy_Data %>%
  separate(Chemical, into = c("Chemical", "Date"), sep = "_") %>% #Seperating chemical and date variable
  mutate(Date = NULL) %>% #Getting rid of date variable
  mutate(Animal = rep(1:54, times = length(folderNames)))
#Ignore warning message... it shous up because separate finds 3 different chunks (because year, month, date are also separated by an undercore... We are discarding the variable anyway)

rm(Baseline_avg, h24_avg, Delta, Delta_1, Delta_2, AllDoses, AllDoses_2)

write_csv(x = Tidy_Data, file = "Data/Alamar_Blue_Tidy_Data_26chems.csv")

####Analysis Bit####

#Outliers
outlierTest <- outlierTest(lm(Fluorescence ~ as.factor(`Dose(mg/L)`) + as.factor(Group), data = Tidy_Data))
#Outlier test object identified these outliers
outliers <- Tidy_Data[c(1268, 771, 799, 115), "Fluorescence"]
Tidy_Data <- replace_with_na(Tidy_Data, outliers)
rm(outliers, outlierTest)

#Mean and StDev of Delta of each Dose
DoseSummary <- Tidy_Data %>%
  group_by(Chemical, `Dose(mg/L)`) %>%
  summarise(
    StDev = sd(Fluorescence, na.rm = TRUE),
    Mean = mean(Fluorescence, na.rm = TRUE),
    .groups = "keep"
  )

#Homogeneity of Variance
LeveneResults <- Tidy_Data %>%
  group_by(Chemical) %>%
  summarise(
    leveneTest(Fluorescence, `Dose(mg/L)`)
    )

VarianceCheck <- if (any(LeveneResults$`Pr(>F)` < 0.05, na.rm = TRUE)) {
  print("Not Equal")
  print(which(LeveneResults$`Pr(>F)`[] < 0.05))
 } else {
    print("Equal")
 }
rm(VarianceCheck)

#ANCOVA
ANCOVA <- Tidy_Data %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~aov(Fluorescence ~ as.factor(`Dose(mg/L)`) + as.factor(Group), data = .))) %>%
  select(-data) 

ANCOVACheck <- ANCOVA %>%
  mutate(model_tidy = map(model, tidy)) %>%
  unnest(model_tidy)

VarianceCheck_2 <- 
  if (any(ANCOVACheck$p.value < 0.05, na.rm = TRUE)) {
  print("Significant")
  print(ANCOVACheck$Chemical[c(which(ANCOVACheck$p.value[] < 0.05))])
} else {
  print("None Significant")
}
rm(ANCOVACheck, VarianceCheck_2)

#ANOVA
ANOVA <- list()
for (i in 1:length(folderNames)) {
  ANOVA[[i]] <- Anova(ANCOVA$model[[i]], type = "3")
}
names(ANOVA) <- folderNames

VarianceCheck_3 <- #Not done
  
####Plots and Visuals####
#Group Variance
Tidy_Data %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, Fluorescence) %>%
  ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Fluorescence)) +
  geom_smooth(aes(group = Group, color = Group), ) +
  geom_point(aes(color = Group)) +
  facet_wrap(~Chemical, scales = "free") +
  xlab("Dose(mg/L)")

#Individual variance (w/ 75% Interquartile Range Removed)
#Replace 75% IQR w/ NA
out <- boxplot.stats(Tidy_Data$Fluorescence)$out
out_index <- which(Tidy_Data$Fluorescence %in% c(out))
replacewithna <- Tidy_Data[out_index,5]
BoxPlotData <- replace_with_na(Tidy_Data, replacewithna)

BoxPlotData %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, Fluorescence) %>%
  ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Fluorescence)) +
  geom_boxplot() +
  facet_wrap(~Chemical, scales = "free") +
  xlab("Dose(mg/L)")

#Dose Variance 
#Columns = mean of all values in dose group 
#points = raw values coloured by replicate group w/ 75% IQR outliers removed
BoxPlotData %>%
  group_by(Chemical, `Dose(mg/L)`) %>%
  summarise(SD = sd(Fluorescence, na.rm = TRUE), Fluorescence = mean(Fluorescence, na.rm = TRUE)) %>%
  ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Fluorescence)) +
  geom_col() +
  facet_wrap(~Chemical, scales = "free") +
  xlab("Dose(mg/L)") +
  geom_point(data = 
               BoxPlotData %>%
               group_by(Chemical) %>%
               select(Chemical, Group, `Dose(mg/L)`, Fluorescence), 
             aes(x = as.factor(`Dose(mg/L)`), 
                 y = Fluorescence, 
                 colour = Group), 
             na.rm = TRUE) +
  geom_errorbar(aes(ymin = Fluorescence - SD, ymax = Fluorescence + SD))