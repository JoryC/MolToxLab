####Libraries####
library(plyr)
library(data.table)
library(abind)
#library(dplyr)
#library(tidyr)
#library(readr)
library(car)
#library(purrr)
library(broom)
library(naniar)
#library(ggplot2)
library(DescTools)
library(tidyverse)

#### Options ####
options(scipen = 9)

####Directory####
getwd() #Output should be "/*/*/MoltToxLab/Alamar_Blue" or something similar
folderNames <-
  list.files("Data/Sub/")#Character list of all folders in /Data/Sub/
folder <- folderNames
files_in_folder <-
  list.files(paste0("Data/Sub/", folder)) #Character list of all files in each folder (alphabetical order)

####Import the Data####


#Easy way to import the data is to just read my final .csv tables...
#all the data in the directory
Tidy_Data <-
  read_csv(file = "Data/Alamar_Blue_Tidy_Data_26chems.csv")
#or
Tylers_Data <- 
  read_csv(file = "Data/Alamar_Blue_Tylers_Data.csv")




#How I produced the final .csv tables...
list <- list() #Create an empty list
for (k in 1:length(files_in_folder)) {
  list[[k]] <-
    fread(
      file = paste0("Data/All/", files_in_folder[[k]]),
      header = FALSE,
      skip = 13,
      nrows = 8,
      select = c(3:11)
    )
}
names(list) <-
  files_in_folder #Give each dataframe/datatable its corresponding name
#Create a quick functioin to remove the pesky empty row
cleanup <- function(x) {
  x <- x[-7, ]
  return(x)
}
#Remove those pesky rows
list <- lapply(list, FUN = cleanup)



#Subset data
#Data is listed in alphabetical order in 'list' object
#REVIEW: Not ideal way to subset data... naming rows by Chemical - Dose can fix this and make subsetting more streamlined
Baseline_1 <-
  c(seq(from = 3, to = 104, by = 4)) #WARNING: seq() to argument is dynamic... will change with more data
Baseline_2 <- c(seq(from = 4, to = 104, by = 4))
h24_1 <-
  c(seq(from = 1, to = 104, by = 4)) #leading h because object can't start with numeric
h24_2 <- c(seq(from = 2, to = 104, by = 4))
#These objects represent the numnber in 1:4 that corresponds to the 24h and baseline plates... so we have 1-4 replicated 26 times (because we have 26 chemicals in the data frame)... 4*26=104
#Create your 4 dataframes to prepare for creating average data sets
list[Baseline_1] #Check and see that the output is Baseline_1.txt files
Base_1_dfl <-
  rbindlist(list[Baseline_1]) #Here we are subsetting the list object which contains data from each plate (replicated twice to control for instrument reading errors). Each subset contains data for different replicates (2) and two time points (baseline and 24h) = 4 different data frames
list[Baseline_2]
Base_2_dfl <- rbindlist(list[Baseline_2])
list[h24_1]
h24_1_dfl <- rbindlist(list[h24_1])
list[h24_2]
h24_2_dfl <- rbindlist(list[h24_2])
#Simple average of the two data sets (Day 0 - Baseline, and Day 1 - 24h)
Baseline_avg <- as.data.frame(rowMeans(abind(Base_1_dfl, Base_2_dfl, along = 3), #Along '3rd dimension'... not row or column
                                       dims = 2)) #what we are doing here is just taking the average of the two replicates from each plate (26 plates - 26 chemicals -- each read twice for each time point in the flurometer) and making a single data frame
h24_avg <- as.data.frame(rowMeans(abind(h24_1_dfl, h24_2_dfl, along = 3),
                                  dims = 2))
Baseline_avg #Combined data frame for every chemical (including blank wells)
h24_avg #Combine data frame
rm(Base_1_dfl,
   Base_2_dfl,
   h24_1_dfl,
   h24_2_dfl,
   Baseline_1,
   Baseline_2,
   h24_1,
   h24_2) #Cleaning up the work environment
#at this point what we have is the unprocessed single-plate reads for both time points (Day 0 and Day 1)... Next we are going to pull out the blank wells and replace any outliers with NAs



#Subset the control data (First 4 cells by column of every 7th row)
#testing and manually verifying with raw data
# test <- Baseline_avg[c(seq(from = 7, to = 182, by = 7)),]
# row.names(test) <- folderNames
# test


#Pulling out the blank 'control' wells for normalization later
#Baseline avg blank well
Base_control <-
  Baseline_avg[c(seq(from = 7, to = 182, by = 7)), c(1:4)] #WARNING: seq() to argument is dynamic... will change with more data
#We want the 7th row from every plate and since the data frame is a collection of every single plate into one big data frame... I use a sequence of 7 by 7 (i.e. 7, 14, 21, 28...) repeating all the way to 182... 7*26=182... These are the rows that are going to be pulled out of the data frame
row.names(Base_control) <- folderNames #naming each row
Base_control
#24h avg blank well
h24_control <- h24_avg[c(seq(from = 7, to = 182, by = 7)), c(1:4)]
row.names(h24_control) <- folderNames
h24_control

#Remove the control row of the *_avg dfs
Baseline_avg <-
  Baseline_avg[-c(seq(from = 7, to = 182, by = 7)), ] #WARNING
#Again, same sequence as we saw above
h24_avg <- h24_avg[-c(seq(from = 7, to = 182, by = 7)), ]


#Checking for Outliers
qplot(as.vector(unlist(Baseline_avg[, c(1:9)]))) #Needs no adjustments
qplot(as.vector(unlist(h24_avg[, c(1:9)]))) #Look at the distribution of samples and try to identify any outliers... Seems to be some outliers larger than 200
qplot(as.vector(unlist(h24_avg[, c(1:9)]))) + xlim(c(-10, 150)) # adjusting to get a better sense of the distribution and include negative values in the plot...
#So... if values are less than 5, or larger than 150 replace with NAs
Baseline_avg[Baseline_avg < 5 |
               Baseline_avg > 150] = NA #Replace and values that are below 5 or above 150 with NA
h24_avg[h24_avg < 5 |
          h24_avg > 150] = NA #Replace and values that are below 5 or above 150 with NA
length(Baseline_avg[is.na(Baseline_avg)]) #Number of NAs in baseline data frame
length(h24_avg[is.na(h24_avg)]) #Number of NAs in 24h data frame
which(is.na(Baseline_avg)) #vector id of which cells were NA
which(is.na(h24_avg)) #vector id of which cells were NA
#Note: We have the same vector id in both the baseline and 24h plates for one well... so we expect to see 5 NAs in the Delta data frame coming up next



# #Testing adding Dose and Chemicals to both baseline and 24h dfs
# #Baseline
# Baseline_avg <- Baseline_avg %>%
#   mutate(Dose = rep(c("Dose_1", "Dose_2", "Dose_3", "Dose_4", "Dose_5", "Dose_6"), times = 26),
#          Chemical = rep(folderNames, times = 1, each = 6)) %>%
#   group_by(Chemical) %>%
#   select(Chemical, Dose, everything())
#
# #24H
# h24_avg <- h24_avg %>%
#   mutate(Dose = rep(c("Dose_1", "Dose_2", "Dose_3", "Dose_4", "Dose_5", "Dose_6"), times = 26),
#          Chemical = rep(folderNames, times = 1, each = 6)) %>%
#   group_by(Chemical) %>%
#   select(Chemical, Dose, everything())



#Creating the normalization factor to be applied to the Delta data frame in the next step
Blank_Delta <- h24_control - Base_control
Blank_Delta
Blank_Delta <- rowMeans(Blank_Delta)
Blank_Delta



#Calculate the Delta of the Day 0 vs Day 1
Delta <- h24_avg - Baseline_avg
colnames(Delta) <-
  c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3", "C_1", "C_2", "C_3")
length(Delta[is.na(Delta)]) #This is the amount of NAs we have in the data frame

#Assign dose names to all corresponding chemicals
Delta_1 <- Delta %>%
  mutate(Dose = rep(
    c("Dose_1", "Dose_2", "Dose_3", "Dose_4", "Dose_5", "Dose_6"),
    times = 26
  ),
  Chemical = rep(folderNames, times = 1, each = 6)) %>%
  group_by(Chemical) %>%
  select(Chemical, everything())

#Tidy up the data
Delta_2 <- Delta_1 %>%
  gather(key = "Replicate",
         value = "Fluorescence",
         A_1:C_3,
         factor_key = TRUE)  %>%
  separate(Replicate,
           into = c("Group", "Replicate"),
           convert = TRUE)
Delta_2

#Assign Dose Values
AllDoses <-
  read.csv(file = "Highest_Dose_Aug.csv",
           skip = 1,
           header = TRUE)
AllDoses_2 <- AllDoses %>%
  gather(key = Dose, value = "Dose(mg/L)", Dose_1:Dose_6)

#Final Data Frame
Tidy_Data <- Delta_2 %>%
  inner_join(AllDoses_2) %>%
  arrange(Chemical) %>%
  separate(Chemical, into = c("Chemical", "Date"), sep = "_") %>% #Seperating chemical and date variable
  mutate(Date = NULL, Dose = NULL) %>% #Getting rid of date variable
  mutate(Animal = rep(1:54, times = length(folderNames))) %>%
  select(Chemical, `Dose(mg/L)`, Animal, everything())
#Ignore warning message... it shous up because separate finds 3 different chunks (because year, month, date are also separated by an undercore... We are discarding the variable anyway)



#Normalize
Tidy_Data <- Tidy_Data %>%
  mutate(Norm_factor = rep(Blank_Delta, each = 54, times = 1)) %>% #Creating a new column to subtract from the Fluorescence column
  mutate(Delta_Fluorescence = Fluorescence - Norm_factor) %>% #Creating the new normalized column and giving it a more informative column name
  mutate(Fluorescence = NULL, Norm_factor = NULL) #Getting rid of the now useless variables

#Where delta fluorescence is NA
temp %>%
  mutate(is_NA = is.na(Delta_Fluorescence)) %>%
  filter(is_NA == TRUE)



#Old Not used... replace with custom outlier function
# outlierTest <- outlierTest(lm(Fluorescence ~ as.factor(`Dose(mg/L)`) + as.factor(Group), data = Tidy_Data))
# #Outlier test object identified these outliers
# outliers <- Tidy_Data[c(1268, 771, 799, 115), "Fluorescence"]
#Data Frame to export with removed outliers
# Tidy_Data <- replace_with_na(Tidy_Data, outliers)

#Subsetting Tyler's Data
Tylers_Data <- Tidy_Data %>%
  filter(Chemical %in% c("BPA", "BPAF", "DES", "EE2", "TGSH"))

#Writing the Data to the 'getwd()/Data/' Directory
write_csv(x = Tidy_Data, file = "Data/Alamar_Blue_Tidy_Data_26chems.csv")
write_csv(x = Tylers_Data, file = "Data/Alamar_Blue_Tylers_Data.csv")

rm(
  Baseline_avg,
  h24_avg,
  Delta,
  Delta_1,
  Delta_2,
  AllDoses,
  AllDoses_2
) #Cleaning up the environment

####Analysis Bit####

#Mean and StDev of Delta of each Dose
DoseSummary <- Tidy_Data %>%
  group_by(Chemical, `Dose(mg/L)`) %>%
  summarise(
    StDev = sd(Delta_Fluorescence, na.rm = TRUE),
    Mean = mean(Delta_Fluorescence, na.rm = TRUE),
    .groups = "keep"
  )
DoseSummary

#Homogeneity of Variance
LeveneResults <- Tidy_Data %>%
  group_by(Chemical) %>%
  summarise(leveneTest(Delta_Fluorescence, as.factor(`Dose(mg/L)`)))

# VarianceCheck <-
#   if (any(LeveneResults$`Pr(>F)` < 0.05, na.rm = TRUE)) {
#     print("Not Equal")
#     print(which(LeveneResults$`Pr(>F)`[] < 0.05))
#   } else {
#     print("Equal")
#   }
#LeveneResults[which(LeveneResults$`Pr(>F)`[] < 0.05),]

#Adding an is.significant column to easily parse significant values in a spreadsheet
LeveneResults <- LeveneResults %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = `Pr(>F)` < 0.05,
    true = TRUE,
    false = FALSE
  ))
LeveneResults

write_csv(x = LeveneResults, file = "Output/Levene_Test_Results.csv")
#rm(VarianceCheck)



#ANCOVA
ANCOVA <- Tidy_Data %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~ aov(
    Delta_Fluorescence ~ as.factor(`Dose(mg/L)`) + as.factor(Group),
    data = .
  ))) %>% #Where 'Group' is the dose group replicate... group A B or C for one of 3 petri dishes in the dose group... This tells us if there were human error in making sure each replicate got the same dose
  select(-data)

ANCOVACheck <- ANCOVA %>%
  mutate(model_tidy = map(model, tidy)) %>%
  unnest(model_tidy)

# VarianceCheck_2 <-
#   if (any(ANCOVACheck$p.value < 0.05, na.rm = TRUE)) {
#   print("Significant")
#   print(ANCOVACheck$Chemical[c(which(ANCOVACheck$p.value[] < 0.05))])
# } else {
#   print("None Significant")
# }

ANCOVACheck <- ANCOVACheck %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = p.value <= 0.05,
    true = TRUE,
    false = FALSE
  ))
ANCOVACheck #Final Data Frame for the ANCOVA... will write later

ANCOVA_Sig_Results <-
  ANCOVACheck[which(ANCOVACheck$p.value[] <= 0.05), ]
ANCOVA_Sig_Results #Just the significant results of the ANCOVA

#ANOVA
ANOVA_Data <- Tidy_Data %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~ aov(
    Delta_Fluorescence ~ as.factor(`Dose(mg/L)`), data = .
  )))

ANOVA_Check <- ANOVA_Data %>%
  mutate(model_tidy = map(model, tidy)) %>%
  unnest(model_tidy)

# VarianceCheck_3 <-
#   if (any(ANOVA_Check$p.value < 0.05, na.rm = TRUE)) {
#     print("Significant")
#     print(ANOVA_Check$Chemical[c(which(ANOVA_Check$p.value[] < 0.05))])
#   } else {
#     print("None Significant")
#   }

ANOVA_Check <- ANOVA_Check %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = p.value < 0.05,
    true = TRUE,
    false = FALSE
  ))
ANOVA_Check #Final Data frame of the ANOVA Results

ANOVA_Sig_Results <-
  ANOVA_Check[which(ANOVA_Check$p.value[] <= 0.05), ]
ANOVA_Sig_Results


#Saving the ANCOVA and ANOVA Results
write_csv(x = ANCOVACheck, file = "Output/ANCOVA_Results.csv") #ANCOVA
write_csv(x = ANOVA_Check, file = "Output/ANOVA_Results.csv") #ANOVA


#PostHoc tests

#Dunnett's Test
Dunnett_results <- Tidy_Data %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~ DunnettTest(
    x = .$Delta_Fluorescence, g = .$`Dose(mg/L)`
  ), data = .)) #Performing the Dunnett's test and saving it is a variable

#Creating list of summaries
#Since the PostHocTest object cannot be coerced to a tidy tibble using tidy()... we got creative
Dunnett_list <-
  list() #What we are trying to do is index the results and see what the significant results were... so we are using a list which can be later coerced into a tibble to easily index...
for (i in 1:length(unique(Tidy_Data$Chemical))) {
  Dunnett_list[[Dunnett_results$Chemical[i]]] <-
    Dunnett_results$model[[i]][["0"]] %>% #Take Dunnett's test results without any of the fancy summary information and shove it into a named list
    as.data.frame() %>% #Coerce to a data frame temporarily so what we can take the row names of the reults and turn them into a variable with rownames_to_column
    rownames_to_column(var = "dose")
}
Dunnett_comb <-
  ldply(Dunnett_list) #this function combines all of the lists together and gives them a variable name according to the chemical
Dunnett_comb$dose = substr(Dunnett_comb$dose,
                           start = 1,
                           stop = nchar(Dunnett_comb$dose) - 2) %>%
  as.numeric() #Here we are fixing the dose column... the dose column has the test dose related to the control... but we just want to see what the test dose is without it giving us redundant information about the comparison to the control for every observation...
Dunnett_comb <- as_tibble(Dunnett_comb) #Coerce to a tidy tibble
Dunnett_comb #Great, a nice tibble that we can export

#Now just to add one more column
Dunnett_comb <- Dunnett_comb %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = pval < 0.05,
    true = TRUE,
    false = FALSE
  ))
Dunnett_comb

write_csv(Dunnett_comb, file = "Output/Dunnett_test_results.csv")

#indexing what the significant results were...
Dunnett_Sig_Results <-
  Dunnett_comb[which(Dunnett_comb$pval <= 0.05), ]
Dunnett_Sig_Results
#Cool!


#Bonferronni Adjustment
#Jasoooooooooooooonnnnnnnnnnn... help



####Plots and Visuals####
#Group Variance
Tidy_Data %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, Delta_Fluorescence) %>%
  ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Delta_Fluorescence)) +
  geom_smooth(aes(group = Group, color = Group), se = FALSE) +
  geom_point(aes(color = Group)) +
  facet_wrap( ~ Chemical, scales = "free") +
  xlab("Dose(mg/L)")

#Playing with geom_smooth
Tidy_Data %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, Delta_Fluorescence) %>%
  ggplot(aes(x = log10(`Dose(mg/L)`), y = Delta_Fluorescence)) +
  geom_smooth(se = FALSE) +
  geom_jitter(width = 0.2) +
  facet_wrap( ~ Chemical, scales = "free") +
  xlab("Dose(mg/L)")

#Individual variance (w/ 75% Interquartile Range Removed)
# #Replace 75% IQR w/ NA
# out <- boxplot.stats(Tidy_Data$Delta_Fluorescence)$out
# out_index <- which(Tidy_Data$Delta_Fluorescence %in% c(out))
# replacewithna <- Tidy_Data[out_index, 5]
# BoxPlotData <- replace_with_na(Tidy_Data, replacewithna)

#Boxplot - Dose effects - Tyler's Data
g <- Tylers_Data %>%
  mutate(Dose = `Dose(mg/L)` * 1000) %>%
  filter(Delta_Fluorescence < 50 & Delta_Fluorescence >= 0) %>%
  group_by(Chemical) %>%
  select(Chemical, Group, Dose, Delta_Fluorescence) %>%
  ggplot(aes(x = as.factor(Dose), y = Delta_Fluorescence)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, height = 0) +
  #geom_smooth(se = FALSE) +
  facet_wrap(~ Chemical, scales = "free", ncol = 2) +
  ylab("Change in Fluorescence") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, hjust = 0)) +
  xlab("Dose (µg/L)")

ggsave(
  units = "in",
  height = 8,
  width = 12,
  filename = "Output/Alamar_Blue_Tyler.pdf",
  plot = g,
  path = getwd(),
  device = "pdf"
)
# ggsave(
#   units = "",
#   height = 800,
#   width = 1200,
#   filename = "Output/Alamar_Blue_Tyler.png",
#   plot = g,
#   path = getwd(),
#   device = "png"
# )

#Dose Variance
#Columns = mean of all values in dose group
#points = raw values coloured by replicate group w/ 75% IQR outliers removed
Tidy_Data %>%
  group_by(Chemical, `Dose(mg/L)`) %>%
  summarise(
    SD = sd(Delta_Fluorescence, na.rm = TRUE),
    Delta_Fluorescence = mean(Delta_Fluorescence, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Delta_Fluorescence)) +
  geom_col() +
  facet_wrap( ~ Chemical, scales = "free") +
  xlab("Dose(mg/L)") +
  geom_jitter(
    data =
      Tidy_Data %>%
      group_by(Chemical) %>%
      select(Chemical, Group, `Dose(mg/L)`, Delta_Fluorescence),
    width = 0.15,
    aes(
      x = as.factor(`Dose(mg/L)`),
      y = Delta_Fluorescence
    ),
    na.rm = TRUE
  ) +
  geom_errorbar(aes(ymin = Delta_Fluorescence - SD, ymax = Delta_Fluorescence + SD))
