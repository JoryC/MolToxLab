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
library(outliers)
library(future) 
library(future.apply) # from original future author
library(furrr) # use future.apply but similar to purrr package
future::plan(multisession, workers = availableCores()) # windows, Mac needs multicore
library(ggrepel)
library(scales)
library(emmeans)
library(multcomp)
library(multcompView)

#### Options ####
options(scipen = 9)

####Directory####
getwd() #Output should be "/*/*/MoltToxLab/Alamar_Blue" or something similar
folderNames <-
  list.files("Data/Sub/")#Character list of all folders in /Data/Sub/
folder <- folderNames
files_in_folder <-
  list.files(paste0("Data/Sub/", folder), full.names = TRUE, pattern = "*.txt") #Character list of all files in each folder (alphabetical order)
fileNames <-   list.files(paste0("Data/Sub/", folder), full.names = FALSE, pattern = "*.txt") #Character list of all files in each folder (alphabetical order)

####Import the Data####


# #Easy way to import the data is to just read my final .csv tables...
# #all the data in the directory
Tidy_Data <-
  read_csv(file = "Data/Alamar_Blue_Tidy_Data_29chems.csv")
Tidy_Data_ol_rm <-
  read_csv(file = "Data/Alamar_Blue_Tidy_Data_29chems_outliers_rm.csv")
# #or
# Tylers_Data <- 
#   read_csv(file = "Data/Alamar_Blue_Tylers_Data.csv")




#How I produced the final .csv tables...



#Read in Data
list <- list() #Create an empty list
for (k in files_in_folder) {
  list[[k]] <-
    fread(
      file = k,
      skip = 12, #We want to just read in the second data bloack fromt he raw .txt friles from the spectrophotometer read-out
      nrows = 8, #Stop reading after the data frame has been read in
      select = c(3:11), #Select just the 9 columns that contain the data and no white space
      header = FALSE,
      data.table = TRUE
    ) # Using fread because it is really quick and multicore. It's faster than read_tsv for example
}
names(list) <- str_split(string = names(list), pattern = "/", simplify = TRUE)[,4]
#Create a quick functioin to remove the pesky empty row between the data and the blank well controls
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
  c(seq(
    from = 3,
    to = length(files_in_folder),
    by = 4
  )) #WARNING: seq() to argument is dynamic... will change with more data
Baseline_2 <- c(seq(
  from = 4,
  to = length(files_in_folder),
  by = 4
))
h24_1 <-
  c(seq(
    from = 1,
    to = length(files_in_folder),
    by = 4
  )) #leading h because object can't start with numeric
h24_2 <- c(seq(
  from = 2,
  to = length(files_in_folder),
  by = 4
))
#These objects represent the numnber in 1:4 that corresponds to the 24h and baseline plates... so we have 1-4 replicated 26 times (because we have 26 chemicals in the data frame)... 4*26=104
#Create your 4 dataframes to prepare for creating average data sets
list[Baseline_1] #Check and see that the output is Baseline_1.txt files
Base_1_dfl <- plyr::ldply(list[Baseline_1]) #Here we are subsetting the list object which contains data from each plate (replicated twice to control for instrument reading errors). Each subset contains data for different replicates (2) and two time points (baseline and 24h) = 4 different data frames
list[Baseline_2]
Base_2_dfl <- plyr::ldply(list[Baseline_2])
list[h24_1]
h24_1_dfl <- plyr::ldply(list[h24_1])
list[h24_2]
h24_2_dfl <- plyr::ldply(list[h24_2])
#Simple average of the two data sets (Day 0 - Baseline, and Day 1 - 24h)
Baseline_avg <-
  rbind(Base_1_dfl, Base_2_dfl) %>% #Change
  mutate(
    Chemical = str_split(
      string = .id,
      pattern = "_",
      simplify = TRUE
    )[, 1],
    Measurement_1 = paste0(
      str_split(
        string = .id,
        pattern = "_",
        simplify = TRUE
      )[, 2],
      "_",
      str_split(
        string = .id,
        pattern = "_",
        simplify = TRUE
      )[, 3]
    )
  ) %>%
  mutate(Measurement = str_remove(string = Measurement_1, pattern = ".txt")) %>%
  mutate(Dose = rep(
    c(
      "Dose_1",
      "Dose_2",
      "Dose_3",
      "Dose_4",
      "Dose_5",
      "Dose_6",
      "Blank"
    ),
    times = length(folderNames) * 2
  )) %>%
  pivot_longer(cols = V3:V11, values_to = "Fluorescence") %>%
  mutate(is_empty = if_else(
    condition = (name %in% c("V7", "V8", "V9", "V10", "V11")) &
      Dose == "Blank",
    true = TRUE,
    false = FALSE
  )) %>%
  filter(is_empty != TRUE) %>%
  select(Chemical, Dose, name, Fluorescence, Measurement) %>%
  pivot_wider(names_from = Measurement, values_from = Fluorescence) %>%
  mutate(Baseline_Fluorescence = (Baseline_1 + Baseline_2) / 2) %>% #Change
  select(-name)

h24_avg <-
  rbind(h24_1_dfl, h24_2_dfl) %>% #Change
  mutate(
    Chemical = str_split(
      string = .id,
      pattern = "_",
      simplify = TRUE
    )[, 1],
    Measurement_1 = paste0(
      str_split(
        string = .id,
        pattern = "_",
        simplify = TRUE
      )[, 2],
      "_",
      str_split(
        string = .id,
        pattern = "_",
        simplify = TRUE
      )[, 3]
    )
  ) %>%
  mutate(Measurement = str_remove(string = Measurement_1, pattern = ".txt")) %>%
  mutate(Dose = rep(
    c(
      "Dose_1",
      "Dose_2",
      "Dose_3",
      "Dose_4",
      "Dose_5",
      "Dose_6",
      "Blank"
    ),
    times = length(folderNames) * 2
  )) %>%
  pivot_longer(cols = V3:V11, values_to = "Fluorescence") %>%
  mutate(is_empty = if_else(
    condition = (name %in% c("V7", "V8", "V9", "V10", "V11")) &
      Dose == "Blank",
    true = TRUE,
    false = FALSE
  )) %>%
  filter(is_empty != TRUE) %>%
  select(Chemical, Dose, name, Fluorescence, Measurement) %>%
  pivot_wider(names_from = Measurement, values_from = Fluorescence) %>%
  mutate(h24_Fluorescence = (`24h_1` + `24h_2`) / 2) %>% #Change
  select(-name)


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
Base_control <- Baseline_avg %>%
  filter(Dose == "Blank") %>%
  group_by(Chemical, Dose) %>%
  summarise(Base_Fluorescence = mean(Baseline_Fluorescence))
  
#24h avg blank well
h24_control <- h24_avg %>%
  filter(Dose == "Blank") %>%
  group_by(Chemical, Dose) %>%
  summarise(h24_Fluorescence = mean(h24_Fluorescence))

#Creating the normalization factor to be applied to the Delta data frame in the next step
Blank_Delta <- full_join(x = Base_control, y = h24_control) %>%
  mutate(Fluorescence = h24_Fluorescence-Base_Fluorescence)
Blank_Delta



#Tidy up the data
Delta <- cbind(Baseline_avg, h24_avg) %>%
  subset(., select = which(!duplicated(names(.)))) %>%
  mutate(Delta_Fluorescence = .$h24_Fluorescence-.$Baseline_Fluorescence)

#Assign Dose Values
AllDoses <-
  read.csv(file = "MetaData.csv",
           skip = 1,
           header = TRUE)
AllDoses_2 <- AllDoses %>%
  gather(key = Dose, value = "Dose(mg/L)", Dose_1:Dose_6)

#Final Data Frame
Tidy_Data <- Delta %>%
  inner_join(AllDoses_2) %>%
  arrange(Chemical) %>%
  #separate(Chemical, into = c("Chemical", "Date"), sep = "_") %>% #Seperating chemical and date variable
  #mutate(Date = NULL, Dose = NULL) %>% #Getting rid of date variable
  mutate(Group = rep(c("A", "B", "C"), each = 3, length.out = length(folderNames)*9*6)) %>%
  mutate(Animal = rep(1:54, times = length(folderNames))) %>%
  select(Chemical, `Dose(mg/L)`, Dose, Animal, Group, everything())
#Ignore warning message... it shous up because separate finds 3 different chunks (because year, month, date are also separated by an undercore... We are discarding the variable anyway)



#Normalize
Blank_Delta <- Blank_Delta %>%
  select(Chemical, Fluorescence) %>%
  rename(Norm_factor = Fluorescence)

Tidy_Data <- Tidy_Data %>%
  group_by(Chemical) %>%
  inner_join(Blank_Delta) %>%
  mutate(Norm_Fluorescence = Delta_Fluorescence - Norm_factor) #Creating the new normalized column and giving it a more informative column name
  #mutate(Fluorescence = NULL, Norm_factor = NULL) #Getting rid of the now useless variables

#Where delta fluorescence is NA
Tidy_Data %>%
  mutate(is_NA = is.na(Delta_Fluorescence)) %>%
  filter(is_NA == TRUE)


#Testing for Outliers...
input_Tidy_Data <- Tidy_Data %>%
  group_by(Chemical) %>%
  select(Norm_Fluorescence) %>%
  mutate(Norm_Fluorescence = round(Norm_Fluorescence, digits = 3)) %>% #Bug with grubbs and large decimal places
  # filter(Chemical %in% c("BPA", "BPAF")) %>% #Testing with smaller data set
  nest()

grubbs.flag <- function(x) {
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  while(pv < 0.05) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(Norm_Fluorescence=x,Outlier=(x %in% outliers)))
}

grubbs_results <- input_Tidy_Data %>%
  mutate(grubbs_results = map(data, ~grubbs.flag(x = .x[[1]])))

Tidy_Data_join <- grubbs_results %>%
  select(-data) %>%
  unnest(cols = c(grubbs_results))

Tidy_Data_test <- Tidy_Data %>%
  mutate(Norm_Fluorescence = round(Norm_Fluorescence, digits = 3)) %>%
  inner_join(Tidy_Data_join, by = c("Chemical", "Norm_Fluorescence"))
Tidy_Data_test <- Tidy_Data_test[-which(duplicated(Tidy_Data_test)),] #Delete duplicates introduced by join

test <- Tidy_Data_test %>%
  group_by(Chemical, Dose) %>%
  tally()

#Outliers flagged
Tidy_Data <- Tidy_Data_test %>%
  mutate(Outlier = Outlier.x)

#Outliers removed
Tidy_Data_ol_rm <- Tidy_Data %>%
  mutate(across(
    .cols = c("Norm_Fluorescence"),
    .fns = ~ replace(x = ., list = Outlier, values = NA) 
  ))

#Subsetting Tyler's Data
Tylers_Data <- Tidy_Data %>%
  filter(Chemical %in% c("BPA", "BPAF", "DES", "EE2", "TGSH"))

#Writing the Data to the 'getwd()/Data/' Directory
write_csv(x = Tidy_Data, file = "Data/Alamar_Blue_Tidy_Data_29chems.csv")
write_csv(x = Tidy_Data_ol_rm, file = "Data/Alamar_Blue_Tidy_Data_29chems_outliers_rm.csv")
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
DoseSummary <- Tidy_Data_ol_rm %>%
  group_by(Chemical, `Dose(mg/L)`) %>%
  summarise(
    StDev = sd(Norm_Fluorescence, na.rm = TRUE),
    Mean = mean(Norm_Fluorescence, na.rm = TRUE),
    .groups = "keep"
  )
DoseSummary

#Test of normalization
qqPlot(Tidy_Data_ol_rm$Norm_Fluorescence)
#ggsave(plot = gg_qqplot, filename = "Output/Images/qqplot.png", device = "png", bg = "white")

#Homogeneity of Variance
LeveneResults <- Tidy_Data_ol_rm %>%
  group_by(Chemical) %>%
  summarise(leveneTest(Delta_Fluorescence, as.factor(`Dose(mg/L)`)))
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
ANCOVA <- Tidy_Data_ol_rm %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~ aov(
    Norm_Fluorescence ~ as.factor(`Dose(mg/L)`) + as.factor(Group),
    data = .
  ))) %>% #Where 'Group' is the dose group replicate... group A B or C for one of 3 petri dishes in the dose group... This tells us if there were human error in making sure each replicate got the same dose
  dplyr::select(model)

ANCOVACheck <- ANCOVA %>%
  mutate(model_tidy = map(model, tidy)) %>%
  unnest(model_tidy)

ANCOVACheck <- ANCOVACheck %>%
  mutate(adj_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(is.significant = if_else(
    condition = adj_p.value <= 0.05,
    true = TRUE,
    false = FALSE
  ))
ANCOVACheck #Final Data Frame for the ANCOVA... will write later

ANCOVA_Sig_Results <-
  ANCOVACheck[which(ANCOVACheck$adj_p.value[] <= 0.05), ]
ANCOVA_Sig_Results #Just the significant results of the ANCOVA
#Saving the ANCOVA Results
write_csv(x = ANCOVACheck, file = "Output/ANCOVA_Results.csv") #ANCOVA




#PostHoc tests

#Dunnett's Test
Dunnett_results <- Tidy_Data_ol_rm %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~ DunnettTest(
    x = .$Norm_Fluorescence, g = .$`Dose(mg/L)`
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
  plyr::ldply(Dunnett_list) #this function combines all of the lists together and gives them a variable name according to the chemical
Dunnett_comb$dose = substr(Dunnett_comb$dose,
                           start = 1,
                           stop = nchar(Dunnett_comb$dose) - 2) %>%
  as.numeric() #Here we are fixing the dose column... the dose column has the test dose related to the control... but we just want to see what the test dose is without it giving us redundant information about the comparison to the control for every observation...
Dunnett_comb <- as_tibble(Dunnett_comb) #Coerce to a tidy tibble
Dunnett_comb #Great, a nice tibble that we can export

#Now just to add one more column
Dunnett_comb <- Dunnett_comb %>%
  mutate(adj_p.value = p.adjust(pval, method = "fdr")) %>%
  mutate(is.significant = if_else(
    condition = pval < 0.05,
    true = TRUE,
    false = FALSE
  ))
Dunnett_comb

#Write to a .csv
write_csv(Dunnett_comb, file = "Output/Dunnett_test_results.csv")

#indexing what the significant results were...
Dunnett_Sig_Results <-
  Dunnett_comb[which(Dunnett_comb$pval <= 0.05), ]
Dunnett_Sig_Results
#Cool!








#TukeyHSD's Test on Groups - testing for 'Group Effects'
Tukey_results <- ANCOVA %>%
  mutate(model = map(model, ~ TukeyHSD(x = .x))) #Performing the Dunnett's test and saving it is a variable

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Following guide here https://schmidtpaul.github.io/DSFAIR/compactletterdisplay.html
#contrasts
#option 1
option1 <- Tukey_results %>%
  mutate(data = map(model, ~ broom::tidy(x = .x)))
#option 2
option2 <- ANCOVA %>%
  mutate(option2 = map(model, ~ emmeans::emmeans(object = .x, specs = "Group"))) %>%
  mutate(option2 = map(option2, ~ pairs(x = .x,  adjust = "Tukey")))
#option 3
option3 <- ANCOVA %>%
  mutate(option3 = map(model, ~ multcomp::cld(emmeans::emmeans(object = .x, specs = "Group"), adjust = "Tukey", details = TRUE, Letters = letters, alpha = 0.05)))

#uniform format
option1 <- option1 %>%
  dplyr::select(data) %>%
  unnest(cols = c(data)) %>% #Tidying Output
  filter(term == "as.factor(Group)") %>%
  mutate(adj.p.value = p.adjust(adj.p.value, method = "fdr"))

option2 <- option2 %>%
  dplyr::select(option2) %>%
  mutate(data = map(option2, ~ as_tibble(x = .x))) %>%
  dplyr::select(data) %>%
  unnest(cols = c(data)) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "fdr"))
  
option3_t <- option3 %>%
  dplyr::select(option3) %>%
  mutate(data = map(option3, pluck(.x = "comparisons"))) %>%
  mutate(data2 = map(data, ~ as_tibble(x = .x))) %>%
  dplyr::select(data2) %>%
  unnest(cols = c(data2)) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "fdr"))
option3 <- option3 %>%
  dplyr::select(option3) %>%
  mutate(data = map(option3, pluck(.x = "emmeans"))) %>%
  mutate(data2 = map(data, ~ as_tibble(x = .x))) %>%
  dplyr::select(data2) %>%
  unnest(cols = c(data2))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#TukeyHSD's Test on Groups - testing for 'Group Effects' continued
Tukey_results <- Tukey_results %>% 
  mutate(data = map(model, ~ broom::tidy(x = .x))) %>%
  dplyr::select(-model) %>% 
  unnest(cols = c(data)) %>% #Tidying Output
  mutate(adj.p.value = p.adjust(adj.p.value, method = "fdr"))

Tukey_Sig_results <- Tukey_results %>%
  filter(adj.p.value < 0.05)

Tukey_join_to_Tidy <- Tukey_results %>%
  filter(term == "as.factor(Group)") %>%
  mutate(Group = str_split(string = contrast, pattern = "-", simplify = TRUE)[,1],
         contrast = str_split(string = contrast, pattern = "-", simplify = TRUE)[,2],
         Group_1 = Group,
         contrast_1 = contrast) %>%
  pivot_longer(cols = c(Group, contrast), values_to = "Group", names_to = "X") %>%
  mutate(contrast = if_else(condition = X == "contrast", true = Group_1, false = contrast_1)) %>%
  dplyr::select(Chemical, Group, contrast, term, null.value, estimate, conf.low, conf.high, adj.p.value)

gg_Tukey <- Tukey_join_to_Tidy %>%
  dplyr::select(Chemical, Group, contrast, adj.p.value)


#Write to a .csv
write_csv(Dunnett_comb, file = "Output/Dunnett_test_results.csv")
write_csv(Tukey_results, file = "Output/Tukey_HSD_results.csv")
write_csv(option_3_t)




#Filtering some interesting results
Interesting_Chemicals <- distinct(as_tibble(c(ANCOVA_Sig_Results$Chemical, Tukey_Sig_results$Chemical, Dunnett_Sig_Results)$.id))


####Plots and Visuals####
#Group Variance v1
# p <- Tidy_Data_ol_rm %>%
#   na.omit() %>%
#   group_by(Chemical) %>%
#   mutate(Norm_Fluorescence_show = as.numeric(
#     between(
#       x = Norm_Fluorescence,
#       left = quantile(Norm_Fluorescence, na.rm = TRUE)[2] - 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE),
#       right = quantile(Norm_Fluorescence, na.rm = TRUE)[4] + 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE)
#     )
#   )) %>%
#   mutate(
#     Norm_Fluorescence_ol_rm = if_else(Norm_Fluorescence_show == 1, true = Norm_Fluorescence, false = NULL)
#   ) %>%
#   select(
#     Chemical,
#     Group,
#     `Dose(mg/L)`,
#     Norm_Fluorescence,
#     Norm_Fluorescence_show,
#     Norm_Fluorescence_ol_rm
#   ) %>%
#   ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Norm_Fluorescence)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_smooth(aes(group = Group, color = Group, y = Norm_Fluorescence_ol_rm), se = FALSE) +
#   geom_jitter(
#     aes(color = Group, alpha = Norm_Fluorescence_show),
#     show.legend = FALSE,
#     height = 0,
#     width = 0.1
#   ) +
#   facet_wrap( ~ Chemical, scales = "free") +
#   xlab("Dose(mg/L)")
# 
# ggsave(
#   width = 1920,
#   height = 1080,
#   units = "px",
#   scale = 3,
#   filename = "Output/Images/Alamar_Blue_Group_Eff.png",
#   plot = p,
#   path = getwd(),
#   device = "png"
# )

#Group Variance v2
y_values_4_geom_text <- Tidy_Data_ol_rm %>%
  group_by(Chemical) %>%
  summarise(max_y = max(Norm_Fluorescence, na.rm = TRUE))
gg_Tukey <- gg_Tukey %>%
  inner_join(y_values_4_geom_text)
gg_data <- Tidy_Data_ol_rm %>%
  na.omit() %>%
  group_by(Chemical) %>%
  mutate(Norm_Fluorescence_show = as.numeric(
    between(
      x = Norm_Fluorescence,
      left = quantile(Norm_Fluorescence, na.rm = TRUE)[2] - 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE),
      right = quantile(Norm_Fluorescence, na.rm = TRUE)[4] + 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE)
    )
  )) %>%
  mutate(
    Norm_Fluorescence_ol_rm = if_else(Norm_Fluorescence_show == 1, true = Norm_Fluorescence, false = NULL)
  ) %>%
  select(
    Chemical,
    Group,
    Dose,
    `Dose(mg/L)`,
    Norm_Fluorescence,
    Norm_Fluorescence_show,
    Norm_Fluorescence_ol_rm
  )
model <- ANCOVA
model_means_cld <- option3 %>%
  group_by(Chemical) %>%
  mutate(Group = fct_reorder(Group, emmean)) %>%
  mutate(cont_group = .group)
# gg_data <- gg_data %>%
#   group_by(Chemical) %>%
#   mutate(Group = fct_relevel(Group, levels(model_means_cld$Group)))

p <- ggplot() +
  #y axis
  scale_y_continuous(
    name = "Fluorescence",
    breaks = pretty_breaks(),
    expand = expansion(mult = (c(0.2, 0.2)))
  ) +
  #x axis
  scale_x_discrete(
    name = "Groups"
  ) +
  #layout
  theme_classic() +
  #black data points
  geom_point(
    data = gg_data,
    aes(y = Norm_Fluorescence, x = Group),
    position = position_nudge(x = -0.2)
  ) +
  #black boxplot
  geom_boxplot(
    data = gg_data,
    aes(y = Norm_Fluorescence, x = Group),
    width = 0.05,
    outlier.shape = NA,
    position = position_nudge(x = -0.1)
  ) +
  #red mean value
  geom_point(
    data = model_means_cld,
    aes(y = emmean, x = Group),
    size = 2,
    color = "red"
  ) +
  #red mean errorbar
  geom_errorbar(
    data = model_means_cld,
    aes(ymin = lower.CL, ymax = upper.CL, x = Group),
    width = 0.05,
    color = "red"
  ) +
  #red letters
  geom_text(
    data = model_means_cld,
    aes(
      y = emmean,
      x = Group,
      label = str_trim(cont_group)
    ),
    position = position_nudge(x = 0.1),
    hjust = 0,
    color = "red"
  ) +
  facet_wrap( ~ Chemical, scales = "free") +
  labs(
    caption = str_wrap("Black dots represent raw data. Red dots and error bars represent (estimated marginal) means ± 95% confidence interval per group. Means not sharing any letter are significantly different by the Tukey-test at the 5% level of significance.", width = 70)
  )
  

#NOTE: Results show to remove 
# 2,4-DMP - Group A
# DES - Group A
# Fadrozole - All groups
# Fenitrothione - Group A
  

# p <- Tidy_Data_ol_rm %>%
#   na.omit() %>%
#   group_by(Chemical) %>%
#   mutate(Norm_Fluorescence_show = as.numeric(
#     between(
#       x = Norm_Fluorescence,
#       left = quantile(Norm_Fluorescence, na.rm = TRUE)[2] - 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE),
#       right = quantile(Norm_Fluorescence, na.rm = TRUE)[4] + 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE)
#     )
#   )) %>%
#   mutate(
#     Norm_Fluorescence_ol_rm = if_else(Norm_Fluorescence_show == 1, true = Norm_Fluorescence, false = NULL)
#   ) %>%
#   select(
#     Chemical,
#     Group,
#     Dose,
#     `Dose(mg/L)`,
#     Norm_Fluorescence,
#     Norm_Fluorescence_show,
#     Norm_Fluorescence_ol_rm
#   ) %>%
#   ggplot(aes(x = Group, y = Norm_Fluorescence)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = as.factor(Dose)),
#               height = 0,
#               width = 0.1) +
#   ggrepel::geom_text_repel(
#     data = gg_Tukey,
#     mapping = aes(
#       x = Group,
#       y = max_y + 8,
#       label = if_else(
#         condition = adj.p.value < 0.05,
#         true = contrast,
#         false = ""
#       )
#     ),
#     inherit.aes = FALSE,
#     direction = "x",
#     color = "black",
#     size = 4
#   ) +
#   
#   xlab("Groups") +
#   scale_y_continuous(expand = expansion(mult = (c(0.2, 0.2)))) +
#   theme_classic()

ggsave(
  width = 1920,
  height = 1080,
  units = "px",
  scale = 3,
  filename = "Output/Images/Alamar_Blue_Group_Eff_Tukey.png",
  plot = p,
  path = getwd(),
  device = "png"
)

# #Playing with geom_smooth
# Tidy_Data %>%
#   group_by(Chemical) %>%
#   select(Chemical, Group, `Dose(mg/L)`, Delta_Fluorescence) %>%
#   ggplot(aes(x = log10(`Dose(mg/L)`), y = Delta_Fluorescence)) +
#   geom_smooth(se = FALSE) +
#   geom_jitter(width = 0.2) +
#   facet_wrap( ~ Chemical, scales = "free") +
#   xlab("Dose(mg/L)")
# 
# #Individual variance (w/ 75% Interquartile Range Removed)
# # #Replace 75% IQR w/ NA
# # out <- boxplot.stats(Tidy_Data$Delta_Fluorescence)$out
# # out_index <- which(Tidy_Data$Delta_Fluorescence %in% c(out))
# # replacewithna <- Tidy_Data[out_index, 5]
# # BoxPlotData <- replace_with_na(Tidy_Data, replacewithna)
# 
# #Boxplot - Dose effects - Tyler's Data
# g <- Tylers_Data %>%
#   mutate(Dose = `Dose(mg/L)` * 1000) %>%
#   filter(Delta_Fluorescence < 50 & Delta_Fluorescence >= 0) %>%
#   group_by(Chemical) %>%
#   select(Chemical, Group, Dose, Delta_Fluorescence) %>%
#   ggplot(aes(x = as.factor(Dose), y = Delta_Fluorescence)) +
#   geom_boxplot(outlier.shape = NA, width = 0.5) +
#   geom_jitter(width = 0.2, height = 0) +
#   #geom_smooth(se = FALSE) +
#   facet_wrap(~ Chemical, scales = "free", ncol = 2) +
#   ylab("Change in Fluorescence") +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 14, hjust = 0)) +
#   xlab("Dose (µg/L)")
# 
# ggsave(
#   units = "in",
#   height = 8,
#   width = 12,
#   filename = "Output/Images/Alamar_Blue_Tyler.pdf",
#   plot = g,
#   path = getwd(),
#   device = "pdf"
# )
# # ggsave(
# #   units = "",
# #   height = 800,
# #   width = 1200,
# #   filename = "Output/Images/Alamar_Blue_Tyler.png",
# #   plot = g,
# #   path = getwd(),
# #   device = "png"
# # )
# 
# #Dose Variance
# #Columns = mean of all values in dose group
# #points = raw values coloured by replicate group w/ 75% IQR outliers removed
# Tidy_Data %>%
#   group_by(Chemical, `Dose(mg/L)`) %>%
#   # filter(Chemical %in% as_vector(Interesting_Chemicals)) %>%
#   summarise(
#     SD = sd(Delta_Fluorescence, na.rm = TRUE),
#     Delta_Fluorescence = mean(Delta_Fluorescence, na.rm = TRUE)
#   ) %>%
#   ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Delta_Fluorescence)) +
#   geom_col() +
#   facet_wrap( ~ Chemical, scales = "free") +
#   xlab("Dose(mg/L)") +
#   geom_jitter(
#     data =
#       Tidy_Data %>%
#       group_by(Chemical) %>%
#       select(Chemical, Group, `Dose(mg/L)`, Delta_Fluorescence),
#     width = 0.15,
#     aes(
#       x = as.factor(`Dose(mg/L)`),
#       y = Delta_Fluorescence
#     ),
#     na.rm = TRUE
#   ) +
#   geom_errorbar(aes(ymin = Delta_Fluorescence - SD, ymax = Delta_Fluorescence + SD))
# 
# #Tweaking geom_col() to only show sig results
# Tidy_Data %>%
#   group_by(Chemical, `Dose(mg/L)`) %>%
#   filter(Chemical %in% as_vector(Interesting_Chemicals)) %>%
#   summarise(
#     SD = sd(Delta_Fluorescence, na.rm = TRUE),
#     Mean_Delta_Fluorescence = mean(Delta_Fluorescence, na.rm = TRUE)
#   ) %>%
#   ggplot(aes(x = as.factor(`Dose(mg/L)`), y = Mean_Delta_Fluorescence)) +
#   geom_col(fill = "coral2") +
#   facet_wrap(~ Chemical, scales = "free") +
#   xlab("Dose (mg/L)") +
#   ylab("Change in Fluorescence (d24h)") +
#   labs(title = "Subset of 8/29 Chemicals ") +
#   geom_jitter(
#     data =
#       Tidy_Data %>%
#       filter(Chemical %in% as_vector(Interesting_Chemicals)) %>%
#       select(Chemical, Group, `Dose(mg/L)`, Delta_Fluorescence),
#     mapping = aes(x = as.factor(`Dose(mg/L)`),
#         y = Delta_Fluorescence),
#     na.rm = TRUE,
#     position = position_jitter(width = 0.2, height = 0, seed = 42069)
#   ) +
#   geom_errorbar(aes(ymin = Mean_Delta_Fluorescence - SD, ymax = Mean_Delta_Fluorescence + SD))

###############################################################################

# Remove those groups that are outliers as seen in Alamar_Blue_Group_Eff_Tukey.png
# 24-DMP - Group A
# DES - Group A
# Fadrozole - All--None?
# Fenitrothione - Group A

Tidy_Data_ol_rm <- read_csv(file = "Data/Alamar_Blue_clean_dat.csv")

y_values_4_geom_text <- Tidy_Data_ol_rm %>%
  group_by(Chemical) %>%
  summarise(max_y = max(Norm_Fluorescence, na.rm = TRUE))
gg_Tukey <- gg_Tukey %>%
  inner_join(y_values_4_geom_text)
gg_data <- Tidy_Data_ol_rm %>%
  na.omit() %>%
  group_by(Chemical) %>%
  mutate(Norm_Fluorescence_show = as.numeric(
    between(
      x = Norm_Fluorescence,
      left = quantile(Norm_Fluorescence, na.rm = TRUE)[2] - 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE),
      right = quantile(Norm_Fluorescence, na.rm = TRUE)[4] + 1.5 * IQR(Norm_Fluorescence, na.rm = TRUE)
    )
  )) %>%
  mutate(
    Norm_Fluorescence_ol_rm = if_else(Norm_Fluorescence_show == 1, true = Norm_Fluorescence, false = NULL)
  ) %>%
  dplyr::select(
    Chemical,
    Group,
    Dose,
    `Dose(mg/L)`,
    Norm_Fluorescence,
    Norm_Fluorescence_show,
    Norm_Fluorescence_ol_rm
  )
model <- ANCOVA
model_means_cld <- option3 %>%
  group_by(Chemical) %>%
  mutate(Group = fct_reorder(Group, emmean)) %>%
  mutate(cont_group = .group)
# gg_data <- gg_data %>%
#   group_by(Chemical) %>%
#   mutate(Group = fct_relevel(Group, levels(model_means_cld$Group)))

p <- ggplot() +
  #y axis
  scale_y_continuous(
    name = "Fluorescence",
    breaks = pretty_breaks(),
    expand = expansion(mult = (c(0.2, 0.2)))
  ) +
  #x axis
  scale_x_discrete(
    name = "Groups"
  ) +
  #layout
  theme_classic() +
  #black data points
  geom_point(
    data = gg_data,
    aes(y = Norm_Fluorescence, x = Group),
    position = position_nudge(x = -0.2)
  ) +
  #black boxplot
  geom_boxplot(
    data = gg_data,
    aes(y = Norm_Fluorescence, x = Group),
    width = 0.05,
    outlier.shape = NA,
    position = position_nudge(x = -0.1)
  ) +
  #red mean value
  geom_point(
    data = model_means_cld,
    aes(y = emmean, x = Group),
    size = 2,
    color = "red"
  ) +
  #red mean errorbar
  geom_errorbar(
    data = model_means_cld,
    aes(ymin = lower.CL, ymax = upper.CL, x = Group),
    width = 0.05,
    color = "red"
  ) +
  #red letters
  geom_text(
    data = model_means_cld,
    aes(
      y = emmean,
      x = Group,
      label = str_trim(cont_group)
    ),
    position = position_nudge(x = 0.1),
    hjust = 0,
    color = "red"
  ) +
  facet_wrap( ~ Chemical, scales = "free") +
  labs(
    caption = str_wrap("Black dots represent raw data. Red dots and error bars represent (estimated marginal) means ± 95% confidence interval per group. Means not sharing any letter are significantly different by the Tukey-test at the 5% level of significance.", width = 70)
  )

print(p)

#You can see from the graph now that the outlier groups have been removed