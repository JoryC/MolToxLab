# Jory Curry

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####                 Instructions               #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Your EDA should include:
# importing
# hygiene and tidying
# transformation, as needed
# wrangling steps, as needed
# several visualizations
# In your script here, include comments throughout that explain key points to the user, including what your research question(s) is/are, where your dataset comes from, and what specific actions you are performing in the script and why.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####                 Version Control            #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Machine Info:
# R version 4.1.3 (2022-03-10) -- "One Push-Up"
# Platform: x86_64-pc-linux-gnu (64-bit) 5.11.0-34-generic / Ubuntu 20.04.4
# Desktop: GNOME 3.36.5
# Hardware: CPU - Intel Core i5-9400F 6 core 4.1GHz / RAM - 15924MiB
# Package Versions:
# lazyeval_0.2.2
# car_3.0-12
# tibble_3.1.6
# ggplot2 Package: Version 3.3.5
# dplyr Package: Version 1.0.8
# plyr_1.8.7
# tidyr: Version 1.2.0
# readr Package: Version 2.1.2
# stringr: 1.4.0
# data.table: 1.14.2
# car: 3.0-12
# DescTools_0.99.44
# Rcurvep_1.2.0
# rlang_1.0.2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                 Background                 ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# The data supplied in the zip folder pertains to behavioural data collected from 5 day-old zebrafish embryos using a Viewpoint ZebraBox and Viewpoint ZebraLab software (https://www.viewpoint.fr/en/p/equipment/zebrabox-for-embryos-or-larvae/publications). This is real data that I collected!! It's bound to have some errors
# There are 26 experiments included, each a unique chemical, and each with 6 different dose groups including a vehicle control group
# Chemical exposures are performed in glass petri dishes (3 dishes per dose group -- A,B,C) for 120 hours
# The behavioural data is collected using an infrared camera over a 50-minute period where the first 20 minutes allow the zebrafish embryos to acclimate to their environment, and for the next 30 miunutes there are 5-minute cycles of light and darkness -- 96-well plate format
# There are 9 fish per dose group and 6 dose groups -- 54 fish per plate
# Zebrafish naturally tend to be more active in the dark
# Neurotoxic chemicals may change the swimming behaviours of fish
# The raw data contains many variables that we will explore once we import the data
# You can browse the meta data in the 'Directory & Meta Data' step of this script
# I come from an environmental ecotoxicology laboratory that explores the feasibility of using acute tests to estimate chronic toxicity for aquatic wildlife (e.g., Early behavioural perturbation analyses, early metabolism analyses, and transcriptomic dose-response modelling)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                   Libraries                ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#library(here)
#library(tidyverse) #Something in the tidyverse package is masking a function that I'm using... Just going to comment this out
library(rlang)
library(Rcurvep)
library(DescTools)
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(data.table)
library(car)
library(lazyeval)
# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
source("Functions/cal_auc_simi_endpoints.R")
source("Functions/behavioural_endpoint_calc.R")
options(scipen = 9) #For displaying scientific notation -- set to print out small numbers w/ no scientific notation

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####           Directory & Meta Data             ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# .XLS files have been converted to .csv files and are included in the directory '~/path/to/files/Data/All'
getwd() # Check out the path to the directory... it should be something like 'path/to/folder/Jory_Behaviour_R_Markdown'
fileNames <- list.files("Data/Raw") #Get the name of each .csv file
fileNames
chemicalNames <- str_split(string = fileNames, pattern = ".csv", simplify = TRUE)[,1] #Identify the chemicals included in the files
chemicalNames
metaData <- read.csv(file = "Data/MetaData.csv") #Import the Meta Data that includes information about the data in the folders
as_tibble(metaData) #CAS is the Chemical Abstract Service, MOA is the Mode of Action. This table includes useful information about the exposure concentrations for each chemical dose in mg/L. We'll use this later to create our final data frame
HighDose <- setNames(metaData[,6], chemicalNames) #subsetting metadata
HighDose

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                 Importing                  #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#Import the raw data as a list
list <- list()
for (i in 1:length(fileNames)) {
  list[[i]] <- read_delim(file = paste0("Data/Raw/", fileNames[[i]]),
                          col_names = TRUE,
                          col_types = "cclnninninninnin", #Where c is character, l is logical, n is numeric, i is integer
                          na = c("", 'NA', "NA", "\t NA", "\tNA") #Specify what to consider NA
  )
}
names(list) <- chemicalNames #Name List objects
list 
colnames(list[[1]]) # Here we can see the names of variables in each data frame... 

# ---------------------------------------------------------------------------- #
# animal -- represents individual animals in the experiment
# Treatment -- The Dose group, and internal replicate (petri dish -- A, B or C)
# an -- Unknown - not useful
# start -- start time of observation in seconds
# end -- end time of observation in seconds
# inact -- Inactivity Counts - the number of times the fish went from being active to inactive over the observation time
# inadur -- Inactivity Duration - the duration of time, in seconds, the fish went from being active to inactive over the observation time (1 minute)
# smlct -- Small Activity Counts - the number of times the fish had a small burst of swim activity over the observation time (1 minute)
# smldur -- Small Activity Duration - the duration of the small burst of swim activity over the observation time (1 minute)
# smldist -- Small Activity Distance - The distance travelled during small bursts of activity
# larct -- Large Activity Counts - the number of times the fish had a large burst of swim activity over the observation time (1 minute)
# lardur -- Large Activity Duration - the duration of the large burst of swim activity over the observation time (1 minute)
# lardist -- Large Activity Distance - The distance travelled during large bursts of activity
# emptyct -- Counts that were neither inactive or active (data recording artifact)
# emptydur -- duration of time fish was neither inactive, or active (data recording artifact) - Almost acts like a confidence value. The closer it is to 60, the more unreliable the data are
# There are 16 Variables, not all are useful
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
#make the data a bit more manageable in one big data frame...
# ---------------------------------------------------------------------------- #

raw_data <- as_tibble(rbindlist(list, idcol = "plate_id")) #name the plates by the chemical names
glimpse(raw_data)
class(raw_data)
raw_data

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####               Investigating                #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

### Checking Expectations ###

## Number of observations (number of rows) ##
# ---------------------------------------------------------------------------- #
96*50*26 # We Expect to see 124,800 rows because we have 96 wells, 50 minutes (one minute per row), and 26 different plates
#Great, now let's investigate how many rows we have
# ---------------------------------------------------------------------------- #

nrow(raw_data) #254,018 rows... 
# Odd, we have double the amount of rows we should have... and more! What's happening?

# ---------------------------------------------------------------------------- #
#We expect to see only 50 observations (1 for each minute) per well (96 wells)
50*60 #We expect the 'end' variable to have a highest value of 3000seconds (50 minutes*60 = 3000seconds)
# ---------------------------------------------------------------------------- #

range(raw_data$end) #ranges from 60 seconds to 3002.2 seconds...
max(raw_data$end) # and some plates exceed 50 minutes by just over 2 seconds...

# ---------------------------------------------------------------------------- #
#Which plates exceed 50 minutes...
# ---------------------------------------------------------------------------- #

head(raw_data[which(raw_data$end > 3000),])
unique(raw_data[which(raw_data$end > 3000),]$plate_id) 
length(unique(raw_data[which(raw_data$end > 3000),]$plate_id)) #It looks like 24 of our 26 plates had observations past our assay cutoff of 50 minutes... some exceeding that time by just a few seconds! These plates will have to be trimmed
length(which(raw_data$end > 3000)) # We can see that 4418 observations surpass the 50 minute stop point of the assay... the extra seconds will have to be discarded
length(which(raw_data$end <= 3000)) #249,600 observations
raw_data[raw_data$end <= 3000,] # Subset of data frame without the extra time observations

# ---------------------------------------------------------------------------- #
254018-4418-124800 # Observations - extra observations - expected observations... hmmm equals 124,800...
#Okay so that accounts for 4418 of the extra rows, what about the other extra 124800 observations
# Since it seems that our observations have double the amount we expect, let's check for repeated rows...
# ---------------------------------------------------------------------------- #

any(duplicated(raw_data) == TRUE) # are any rows duplicated exactly?... this won't work so a deeper dive is necessary
# FALSE...
head(raw_data) # We can see that 'animal' is repeated twice, and 'an' contains a FALSE and TRUE for the duplicated 'animal' observation... all of the values seem to be duplicated also... with the exception of 'emptyct' and 'emptydur'... to clean up the duplicates, delete the duplicated observations that are TRUE under the 'an' variable.

# ---------------------------------------------------------------------------- #
# !!! CLEANING !!! 
# Two tasks to do to clean up nrows: 1) Delete observations past 50 mins, 2) Delete duplicated rows #
# ---------------------------------------------------------------------------- #

#Making a temporary data frame to keep working on investigating the data... we'll clean the raw_data more concisely in the next section of the script:
temp <- raw_data[raw_data$end <= 3000,] # Deleting observations past 50 mins (3000 seconds)
temp <- temp %>%
  filter(an == FALSE) %>% # Removing duplicate rows
  select(-c(an)) # These variables are not very informative
nrow(temp) #Good, now we have the expected amount of observations for 26 experiments of a 96-well plate, 50-minute assay (124,800).

# ---------------------------------------------------------------------------- #
# Number of observations per treatment group #
# ---------------------------------------------------------------------------- #

temp %>%
  na.omit() %>% #NAs represent wells with no fish in them (empty wells)
  group_by(plate_id, Treatment) %>% # Each treatment has 3 fish per treatment, and 50 observations (50 minutes)
  tally() %>% # So we expect to see every treatment with 150 observations total
  summarise(mean_n = mean(n)) #if we take the mean of all treatments, they should all theoretically be exactly 150
#Looks good for every plate... but those NAs are a problem...


## NAs ##
# ---------------------------------------------------------------------------- #
# !!! CLEANING !!!
# We expect wells 10:12, 22:24, 34:36, 46:48, 58:60, 70:96 to all be NA in the 'Treatment' Column because these are all empty wells
(3*6) + (12*2) #Each plate has 42 empty wells
42*50*26 #This means that all NA treatments should be 54,600 observations long
# ---------------------------------------------------------------------------- #

naTemp <- temp[which(is.na(temp$Treatment)),]
nrow(naTemp) #Good, it's 54,600

nrow(temp) - nrow(naTemp) # 70,200... observations, minus empty wells

temp <- temp[-which(is.na(temp$Treatment)),] #So 124,800 observations - 54,600 NA observations = 70,200
nrow(temp) #Good



## Suspicious Values ##
# ---------------------------------------------------------------------------- #
# !!! CLEANING !!!
#If the 'emptydur' variable is maxed out at 60seconds, that means that for that entire one minute, the software was not able to detect a fish in the well at all with its infrared camera. However, we know that those wells have fish in them... they could have been active or inactive during this technical glitch. Therefore, turning all the values in to NAs in this case makes sense because we can not be confident that they are true zeroes in any of the variables.
# ---------------------------------------------------------------------------- #

#Empty Duration
#Determining a cutoff value for sketchy observations
qplot(data = temp, x = emptydur) # We can see that some wells are being recorded as empty with no fish... let's see what fish... It look like anything above 3 seconds of empty duration is an error from this figure
sketchyCutoff <- (60*0.05) #My arbitrary cutoff value
temp[which(temp$emptydur >= sketchyCutoff),] #Ah yes, suspicious fish...  this could be due to the fish being dead or artifacts of the recording software... we can see from this tibble that we have 1,641 rows

#Exploring the Sketchy values
nSketchy <- nrow(temp[which(temp$emptydur >= sketchyCutoff),]) 
nSketchy #Number of sketchy observations... 1,641
qplot(data = temp[-which(temp$emptydur >= sketchyCutoff),], x = emptydur) # It appears that a cutoff of 3 (60*0.05) is good

#Removing the sketchy values
temp <- temp %>%
  mutate(Suspicious = if_else(condition = emptydur >= sketchyCutoff, true = TRUE, false = FALSE)) #Flagging any suspicious false zeroes... I used 3 empty seconds (5% of 60) as my arbitrary cutoff
temp[temp$Suspicious, c("inact", "inadur", "inadist", "smlct", "smldur", "smldist", "larct", "lardur", "lardist")] = NA #Turn all of the rows that are suspicious inta NAs but only in the listed variables
temp %>%
  filter(Suspicious == TRUE) #There, now we can see that those sketchy zero observations have been turned into NAs for each of the variables we wanted. 1,641 observations turned into NAs!

# ---------------------------------------------------------------------------- #
#Overall, the animal recording set up had an approximate failure-rate of 2% -- This is the approximate percentage of time the infrared camera failed to detect an animal when it was present
paste0(round((nSketchy/nrow(temp))*100), "%", " ", "Failure-rate") 
# ---------------------------------------------------------------------------- #


## Variable Observations (makes sense?), Distributions and Outliers ##
# ---------------------------------------------------------------------------- #
#Next up is to check that each variable has a range of values that makes sense. 
# ---------------------------------------------------------------------------- #
#26 chemicals/plates
length(unique(temp$plate_id)) #Good

#There are 96 wells in the plate...
unique(temp$animal)
length(unique(temp$animal)) #Okay we have 54 'animals' (9 fish per dose group, 6 dose groups = 54 fish), good

#6 Doses per chemical X 3 Replicates (A,B,C) per treatment = 18 'Treatments'
unique(temp$Treatment) #NAs are still included
length(unique(na.omit(temp$Treatment))) #Omit the NAs and count 'Treatment' groups... 18 good

# ---------------------------------------------------------------------------- #
# Time #
# ---------------------------------------------------------------------------- #

range(temp$start) #Should be 0 to 2,940 seconds... good
range(temp$end) #60 - 3,000 seconds

# ---------------------------------------------------------------------------- #
# Counts #
# ---------------------------------------------------------------------------- #

# -- Inactive
range(na.omit(temp$inact)) #Inactivity -- This number can vary quite a bit but they should all be integers
is.integer(na.omit(temp$inact)) #Good - TRUE
# - Inactive Counts Quick Plot
qplot(data = temp, x = inact) # Normal Dist

# -- Small activity
range(na.omit(temp$smlct)) #Small Activity -- This number can vary quite a bit
is.integer(na.omit(temp$smlct)) #Good - TRUE
# - Small Counts Quick Plot
qplot(data = temp, x = smlct) # Looks Normal - Maybe a couple small outliers but nothing crazy

# -- Large activity
range(na.omit(temp$larct)) #Large Activity -- This number can vary quite a bit but they should all be integers
is.integer(na.omit(temp$larct)) #Good - TRUE
# - Large Counts Quick Plot
qplot(data = temp, x = larct) # Counts seem to cluster around 0 counts or just above 100 counts


# ---------------------------------------------------------------------------- #
# Duration #
# !!! CLEANING !!!
# This variable should not exceed 60 seconds for any variable, and when we sum all the variables up together from the same row, they should equal to 60 seconds total
# ---------------------------------------------------------------------------- #

# -- Inactivity
range(na.omit(temp$inadur)) #Good
# - Inactive Duration Quick Plot
qplot(data = temp, x = inadur) # Normal distribution but a bit skewed

# -- Small Activity
range(na.omit(temp$smldur)) #Good
# - Small Duration Quick Plot
qplot(data = temp, x = smldur) # Looks Normal

# -- Large Activity
range(na.omit(temp$lardur)) #Good
# - Large Duration Quick Plot
qplot(data = temp, x = lardur) # large duration again seems to cluster closer to 0 or ~15. Outliers past 35 seconds
length(which(temp$lardur >= 46)) #5 outliers - could use a bit of cleaning potentially
temp[which(temp$lardur >= 46),] #Seems they are all from the same animal, plate and dose... We'll have to filter the data to that plate to see the distribution of that plate to see if it truly might be an outlier... This might be useful data (especially considering that it is in the highest dose group) and we don't want to just carelessly get rid of it...
temp %>%
  filter(plate_id == "BPA") %>% #Subetting to just BPA...
  qplot(data = ., x = lardur) # Although it does seem like an outlier, we can't know for sure if the observations are observable effects from the treatment or just random noise
temp %>%
  filter(plate_id == "BPA") %>%
  mutate(Treatment = rep(
    c(
      "Dose1_A","Dose1_B","Dose1_C","Dose2_A","Dose2_B","Dose2_C",
      "Dose3_A","Dose3_B","Dose3_C","Dose4_A","Dose4_B","Dose4_C",
      "Dose5_A","Dose5_B","Dose5_C","Control_A","Control_B","Control_C"
    ),
    times = (nrow(.)/(18)/3), #repeat 2,700/18 = 150 times (number of observations divided by number of treatments = number of treatment observations) then divide that by 3 because each will be repeating 3 times
    each = 3
  )) %>% #Separating Dose from Group letter with an underscore
  separate(col = Treatment, # Column we want to separate
           into = c('Dose', 'Group'), #Into these two column names
           convert = TRUE) %>%
  ggplot(aes(x = lardur)) +
  geom_histogram() +
  facet_wrap(~Dose) # Plot all the distributions of this variable by dose and we will see how the distributions shift from dose to dose
#We see that the highest dose (Dose 1) has more higher value outliers than the rest of the doses and less zero value large activity duration observations... this leads me to believe that this might be an effect of the treatment on behaviour
#Verdict is: Keep the 'outliers' for large activity duration observations

# -- Total Activity Duration
temp <- temp %>%
  mutate(totaldur = emptydur + inadur + smldur + lardur) # Creating another column to test and see if sum of all durations equals 60 for each row
range(na.omit(temp$totaldur)) #Looks okay, but at least one observation has an extra 0.2 second... We can also see that some values are just below at 59.8 seconds... This could just just be an artifact of the recording software. But it is promising that all of these values add up to approximately 60 seconds which tells us that the animal-tracking software is working most of the time except for those instances where we had 'suspicious zeroes'...
qplot(data = temp, x = totaldur) # What we expected to see

# ---------------------------------------------------------------------------- #
# Distance #
# !!! CLEANING !!!
#This variable can also vary quite a bit...
# ---------------------------------------------------------------------------- #

# -- Inactivity
range(na.omit(temp$inadist)) #Kind of a dumb variable here... If a fish is inactive, it shouldn't be moving anywhere at all. Perhaps the reason we see a range of 0-6.1 is because the distance of the activity was below a certain threshold and so was scored as an inactive distance... either way, this variable can go
# - Inactive Distance Quick Plot
qplot(data = temp, x = inadist) # All 0, except some a few strange outliers
length(which(temp$inadist > 0)) #59 values outside of the expected range
qplot(data = temp[-which(temp$inadist > 0),], x = inadist) # Beautiful gray rectangle... that's what we want to see
# This variable will be removed

# -- Small Activity
range(na.omit(temp$smldist)) #Okay
# - Small activity Distance Quick Plot
qplot(data = temp, x = smldist) #Normal dist. with a few outliers it seems past 1250
length(which(temp$smldist >= 1250)) #15 outliers - nothing major
qplot(data = temp[-which(temp$smldist >= 1250),], x = smldist) # Still a nice normal distribution just looks a bit better
# Exploring those 15 outliers
temp %>%
  mutate(Treatment = rep(
    c(
      "Dose1_A","Dose1_B","Dose1_C","Dose2_A","Dose2_B","Dose2_C",
      "Dose3_A","Dose3_B","Dose3_C","Dose4_A","Dose4_B","Dose4_C",
      "Dose5_A","Dose5_B","Dose5_C","Control_A","Control_B","Control_C"
    ),
    times = (nrow(.)/(18)/3), #repeat 70,200/18 = 150 times (number of observations divided by number of treatments = number of treatment observations) then divide that by 3 because each will be repeating 3 times
    each = 3
  )) %>% #Separating Dose from Group letter with an underscore
  separate(col = Treatment,
           into = c('Dose', 'Group'), # New Column names
           convert = TRUE) %>%
  filter(smldist >= 1250) %>%
  group_by(Dose, plate_id) %>%
  tally() #First let's see what dose groups the outliers are from... We can see that 11/15 observations were from the first dose, 3/15 were fromt he thrid dose group and 1/15 were fromt he last dose group...
# Those more potent chemicals seem to be TGSH and EE2... This is consistent with the background information on these chemicals. EE2 is extremely potent, and TGSH has very little toxicological data available at all in the literature (especially for fish -- zero literature)
# This leads me to believe that these 'outliers' may actually be due to the effects of chemical exposures and so, they should be kept in the analysis

# -- Large Activity Distance
range(na.omit(temp$lardist)) #Okay
#However, we expect to see fish moving at least a little bit if they are alive... 
#Large Distance
qplot(data = temp, x = lardist) # Lots of outliers skewing the distribution just above ~2100
length(which(temp$lardist >= 2100)) #158 outliers - could use some cleaning also
qplot(data = temp[-which(temp$lardist >= 2100),], x = lardist) #Looks much better with outliers removed
mean(filter(.data = temp, lardist >= 2100)$inadur) # If the fish were active for a very large distance, we should expect them to be inactive for a short duration. We can see that the mean amount of time an 'outlier' observation was inactive is 4.9 seconds.
mean(filter(.data = temp, lardist <= 2100 & lardist > 0)$inadur) # We can see that for the mean amount of time the other observations spent inactive is 19.2 seconds
#Therefore, the 'outliers' are just hyperactive observations. The cause of the hyperactivity could be a results of the chemical treatment
temp %>%
  filter(lardist >= 2100) %>%
  group_by(plate_id) %>%
  tally() # We see that about 19 chemicals can be found in this subset of hyperactive data... However, it is interesting that BPA, DES, Fadrozole, Flutamide and Vinclozolin are being proportionally more represented that the rest of the subset of data.
# BPA is estrogenic, DES is estrogenic, Vinclozolin is an Androgen receptor antagonist, Flutamide is an Androgen receptor antagonist, and Fadrozole is an Aromatase inhibitor
temp %>%
  mutate(Treatment = rep(
    c(
      "Dose1_A","Dose1_B","Dose1_C","Dose2_A","Dose2_B","Dose2_C",
      "Dose3_A","Dose3_B","Dose3_C","Dose4_A","Dose4_B","Dose4_C",
      "Dose5_A","Dose5_B","Dose5_C","Control_A","Control_B","Control_C"
    ),
    times = (nrow(.)/(18)/3), #repeat 70,200/18 = 150 times (number of observations divided by number of treatments = number of treatment observations) then divide that by 3 because each will be repeating 3 times
    each = 3
  )) %>% #Separating Dose from Group letter with an underscore
  separate(col = Treatment,
           into = c('Dose', 'Group'), # New Column names
           convert = TRUE) %>%
  filter(lardist >= 2100) %>%
  group_by(plate_id, Dose) %>%
  tally() %>%
  print(n = 37)
# We can see from this print out to the console that BPA has its effects at Dose 1 (highest dose), DES has most of its effects in Dose 2, Fadrozole has most of its effects in Dose 1, Flutamide in dose 4, and Vinclozolin mostly in dose 1
# Since there is a discernible pattern, I don't believe these 'outliers' warrant being removed from the analysis

# -- Total Activity Distance
temp <- temp %>%
  mutate(totaldist = smldist + lardist)
range(na.omit(temp$totaldist)) # We see that total dist ranges from 0 to 9067.5... The zero tells us that sometimes the fish don't move at all in the one-minute observation! 
temp %>%
  filter(totaldist == 0) # For those instances where the distance traveled is zero, the inactivity duration or the empty duration explain the lack of distance traveled... they were inactive or could not be detected (likely due to inactivity)
qplot(data = temp, x = totaldist) #Has outliers above 4375
length(which(temp$totaldist >= 4375)) #10 Outliers -- approximately the same as large distance outliers
qplot(data = temp[-which(temp$totaldist >= 4375),], x = totaldist) #Looks much better -- Nice Normal Distribution but with a long tail
temp %>%
  mutate(Treatment = rep(
    c(
      "Dose1_A","Dose1_B","Dose1_C","Dose2_A","Dose2_B","Dose2_C",
      "Dose3_A","Dose3_B","Dose3_C","Dose4_A","Dose4_B","Dose4_C",
      "Dose5_A","Dose5_B","Dose5_C","Control_A","Control_B","Control_C"
    ),
    times = (nrow(.)/(18)/3), #repeat 70,200/18 = 150 times (number of observations divided by number of treatments = number of treatment observations) then divide that by 3 because each will be repeating 3 times
    each = 3
  )) %>% #Separating Dose from Group letter with an underscore
  separate(col = Treatment,
           into = c('Dose', 'Group'), # New Column names
           convert = TRUE) %>%
  filter(totaldist >= 4375) %>%
  group_by(plate_id, Dose) %>%
  tally()
# We see that about 3 chemicals can be found in this subset of hyperactive data...
# BPA is estrogenic, Vinclozolin is an Androgen receptor antagonist, and Fadrozole is an Aromatase inhibitor

#Therefore, after a closer look at the data, the distributions, and potential outliers, it has been determined that sufficient outliers have been removed from the larger data set...

# -----------------------------------------------------------------------------#
#Looking at each set of data individually will reveal any further potential outliers... in total distance variable
# -----------------------------------------------------------------------------#

#Write a function to observe the distributions of different variables in different chemicals
distVis <- function (df, chem, var, time_start = 0) {
df %>%
  filter(plate_id == chem, start >= time_start*60) %>% #Subetting to chemical of choice
  mutate(Treatment = rep(
    c(
      "Dose1_A","Dose1_B","Dose1_C","Dose2_A","Dose2_B","Dose2_C",
      "Dose3_A","Dose3_B","Dose3_C","Dose4_A","Dose4_B","Dose4_C",
      "Dose5_A","Dose5_B","Dose5_C","Control_A","Control_B","Control_C"
    ),
    times = (nrow(.)/(18)/3), #repeat 2,700/18 = 150 times (number of observations divided by number of treatments = number of treatment observations) then divide that by 3 because each will be repeating 3 times
    each = 3
  )) %>% #Separating Dose from Group letter with an underscore
  separate(col = Treatment, # Column we want to separate
           into = c('Dose', 'Group'), #Into these two column names
           convert = TRUE) %>%
  ggplot(aes(x = .data[[var]])) +
  geom_histogram() +
  xlab(var) +
  ggtitle(chem) +
  facet_wrap(~Dose)
}

#Save images to identify outliers visually
for (i in 1:length(chemicalNames)) {
tempDist <- distVis(temp, chem = chemicalNames[i], var = "totaldist")
ggsave(plot = tempDist, path = paste0(getwd(), "/Output", "/", chemicalNames[i]), filename = paste0("Distribution_of_Total_Distance_", chemicalNames[i], ".png"), device = "png", width = 1200, height = 800, units = "px")
}
#Summary of outliers: 24DMP > 3000, BPA > 3000, Clofib > 3000, EE2 > 2500, Fadrozole > 4500, Fluoxetine > 2500, Vinclozolin > 4000

#Save last 30 minute images to identify outliers visually
for (i in 1:length(chemicalNames)) {
  tempDist <- distVis(temp, chem = chemicalNames[i], var = "totaldist", time_start = 20)
  ggsave(plot = tempDist, path = paste0(getwd(), "/Output", "/", chemicalNames[i]), filename = paste0("No_Acclimation_Distribution_of_Total_Distance_", chemicalNames[i], ".png"), device = "png", width = 1200, height = 800, units = "px")
}
#Summary of outliers: TGSH > 2000, Clofib > 2000, 17beta_estradiol > 2500, 24DMP > 2500, DES > 2500, DMF > 2500, EE2 > 2500, Fluoxetine > 2500, BPA > 3000, Fadrozole > 3500, Vinclozolin > 4000

out_1 <- temp %>%
  filter(plate_id == "Vinclozolin") %>%
  filter(totaldist >= 4000)
out_2 <- temp %>%
  filter(plate_id == "Fadrozole") %>%
  filter(totaldist >= 3500)
out_3 <- temp %>%
  filter(plate_id == "BPA") %>%
  filter(totaldist >= 3000)
out_4 <- temp %>%
  filter(plate_id %in% c("Fluoxetine", "EE2", "DMF", "DES", "24DMP", "17beta_estradiol")) %>%
  filter(totaldist >= 2500)
out_5 <- temp %>%
  filter(plate_id %in% c("TGSH", "Clofibric_Acid")) %>%
  filter(totaldist >= 2000)
out_totaldist <- rbind(out_1, out_2, out_3, out_4, out_5)

temp %>%
  mutate(outliers = if_else(condition = (do.call(paste0, .) %in% do.call(paste0, out_totaldist)), true = TRUE, false = FALSE)) %>%
  select(plate_id, Treatment, end, totaldist, outliers) %>%
  View()

# -----------------------------------------------------------------------------#
#Looking at each set of data individually will reveal any further potential outliers... in active duration variable
# -----------------------------------------------------------------------------#

#Creating activedur variable
temp <- temp %>% mutate(activedur = smldur + lardur)

#Save last 30 minute images to identify outliers visually
for (i in 1:length(chemicalNames)) {
  tempDist <- distVis(temp, chem = chemicalNames[i], var = "activedur", time_start = 20)
  ggsave(plot = tempDist, path = paste0(getwd(), "/Output", "/", chemicalNames[i]), filename = paste0("No_Acclimation_Distribution_of_Active_Duration_", chemicalNames[i], ".png"), device = "png", width = 1200, height = 800, units = "px")
}
#Summary of outliers: None

distVis(df = temp, chem = "TGSH", var = "inact", time_start = 20)
distVis(df = temp, chem = "BPAF", var = "inact", time_start = 20)
distVis(df = temp, chem = "DMSO", var = "activedur", time_start = 20)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                  Cleaning                  #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# ---------------------------------------------------------------------------- #
# Adding dose information from the meta data to the data frame
#Tasks to do to clean up data set: 1) Delete observations past 50 mins, 2) Delete duplicated rows, *3) Removing Empty wells (Treatment - NAs)
# *3) # We expect wells 10:12, 22:24, 34:36, 46:48, 58:60, 70:96 to all be NA in the 'Treatment' Column because these are all empty wells 
(3*6) + (12*2) #Each plate has 42 empty wells
42*50*26 #This means that all NA treatments should be 54,600 observations long (124,800 - 54,600 = 70,200)
# The goal is to have 70,200 rows or observations
# ---------------------------------------------------------------------------- #

doseData <- metaData %>%
  select(plate_id, Dose1:Control) %>%
  gather(key = Dose, value = "Dose_mg_L", Dose1:Control)
doseData #Now this is ready to join to the data frame

treatment_groups <- c("Dose1_A","Dose1_B","Dose1_C",
                      "Dose2_A","Dose2_B","Dose2_C",
                      "Dose3_A","Dose3_B","Dose3_C",
                      "Dose4_A","Dose4_B","Dose4_C",
                      "Dose5_A","Dose5_B","Dose5_C",
                      "Control_A","Control_B","Control_C") #Defining all of the dose groups and replicate groups separated by an '_' to make the strings a bit easier to manipulate

fishBehavDat <- raw_data %>% #254,018 Observations
  filter(end <= 3000) %>% #249,600 Observations: 1) Deleting observations past 50 mins (3000 seconds)
  filter(an == FALSE) %>% #124,800 Observations: 2) Removing duplicate rows
  filter(!is.na(Treatment)) %>% #70,200 Observations: *3) Remove the empty wells
  mutate(Treatment = rep( 
    treatment_groups,
    times = (nrow(.) / (18) / 3),
    each = 3 #repeat 70,200/18 = 150 times (number of observations divided by number of treatments = number of treatment observations) then divide that by 3 because each will be repeating 3 times
  )) %>% 
  separate(col = Treatment,
           into = c('Dose', 'Group'), # New Column names
           convert = TRUE) %>% #Separating Dose from Group letter with an underscore
  inner_join(doseData) #Finally, Adding the dose information

# ---------------------------------------------------------------------------- #
# Tasks to do to clean up data set: *4) Mark recording errors and/or dead fish as NA, 5) Create a 'time' variable in minutes instead of seconds, create a variable that represents the total active duration of fish during the 1-minute observation, create a variable a represents the total distance fish traveled during the 1-minute observation, and create another logical variable that tells us if the observation is from the vehichle control group or not, splitting up the Dose from the replicate group (A,B,C), removing redundant 'animal' names, creating a 'replicate' variable that tells us what biological replicate # it is for the dose group, and adding a 'embryo_id' variable for further downstream analysis, 6) Remove variables that will not be used for further analysis and reorder variables
# *4) If the 'emptydur' variable is maxed out at 60seconds, that means that for that entire one minute, the software was not able to detect a fish in the well at all with its infrared camera. However, we know that those wells have fish in them... they could have been active or inactive during this technical glitch. Therefore, turning all the values in to NAs in this case makes sense because we can not be confident that they are true zeroes in any of the variables.
# ---------------------------------------------------------------------------- #

sketchyValueCutoff <- (60*0.05) #My arbitrary cutoff value - 3 seconds of 'empty duration' is all I will tolerate

fishBehavDat <- fishBehavDat %>%
  mutate(sketchy = if_else( #70,200 Observations: *4) Turn all of the rows that are sketchy inta NAs but only in the listed variables
    condition = emptydur >= sketchyValueCutoff | inadist > 0,
    true = TRUE,
    false = FALSE
  )) %>%
  mutate(across(
    .cols = c(
      "inact","inadur","inadist",
      "smlct","smldur","smldist",
      "larct","lardur","lardist"
    ),
    .fns = ~ replace(x = ., list = sketchy, values = NA)
  )) %>% # 5) Adding some additional variables
  mutate(
    time_end = end / 60,
    activedur = smldur + lardur,
    totaldist = smldist + lardist,
    is_VC = if_else(
      condition = Dose == "Control",
      true = 1,
      false = 0),
    animal = as.numeric(str_extract_all(animal, "[0-9]+")),
    replicate = rep(c(1:9), times  = nrow(.)/9),
    embryo_id = paste(Dose, replicate, sep = "_")
  ) %>%
  select(plate_id, embryo_id, replicate, Dose_mg_L, is_VC, Group, time_end, inact, smlct, larct, inadur, smldur, lardur, activedur, smldist, lardist, totaldist) # 6) Changing order of variables and removing those that are not useful

fishBehavDat %>% filter(plate_id == "34DCA") #In the console output, see that Dose 1 of 34DCA was full of dead fish... these have been replaced with NA

# ---------------------------------------------------------------------------- #
# Defining outliers in the total distance variable
# ---------------------------------------------------------------------------- #

#Slow way of subsetting data into outliers
#Outliers were determined visually from the 'No_Acclimation_Distributions' in the 'Output' folder
out_1 <- fishBehavDat %>%
  filter(plate_id == "Vinclozolin") %>%
  filter(totaldist >= 4000)
out_2 <- fishBehavDat %>%
  filter(plate_id == "Fadrozole") %>%
  filter(totaldist >= 3500)
out_3 <- fishBehavDat %>%
  filter(plate_id == "BPA") %>%
  filter(totaldist >= 3000)
out_4 <- fishBehavDat %>%
  filter(plate_id %in% c("Fluoxetine", "EE2", "DMF", "DES", "24DMP", "17beta_estradiol")) %>%
  filter(totaldist >= 2500)
out_5 <- fishBehavDat %>%
  filter(plate_id %in% c("TGSH", "Clofibric_Acid")) %>%
  filter(totaldist >= 2000)
#Bind list of outliers
out_totaldist <- rbind(out_1, out_2, out_3, out_4, out_5)

#Changing the data frame - Turning outliers into NAs
# fishBehavDat <- fishBehavDat %>%
#   mutate(outliers = if_else(condition = (do.call(paste0, .) %in% do.call(paste0, out_totaldist)), true = TRUE, false = FALSE)) %>%
#   mutate(across(
#     .cols = c(
#       "inact","inadur",
#       "smlct","smldur","smldist",
#       "larct","lardur","lardist",
#       "activedur", "totaldist"
#     ),
#     .fns = ~ replace(x = ., list = outliers, values = NA)
#   )) %>%
#   select(-outliers)
# fishBehavDat

# ---------------------------------------------------------------------------- #
#Creating a nested list of data frames (1 chemical per data frame) to use later for an analysis pipeline
# ---------------------------------------------------------------------------- #

#Specify what variable you want to analyze
lfishBehavDat <- fishBehavDat %>%
  mutate(value = totaldist) %>% #We will specify total activity distance
  filter(time_end >= 21) %>% # For the analysis, we don't want the 20 minute acclimation period
  mutate(time_end = time_end-20) %>% # Display time 1-30 minutes
  select(plate_id, embryo_id, is_VC, time_end, value) # Pipeline is picky about the data frames
lfishBehavDat <- split(as.data.frame(lfishBehavDat), ~ plate_id)

# ---------------------------------------------------------------------------- #
# Exporting the clean Behavioural Data
# ---------------------------------------------------------------------------- #

write.csv(x = fishBehavDat, file = paste0(getwd(), "/Data", "/Fish_Behaviour_Data.csv"))
#write.csv(x = fishBehavDat, file = paste0(getwd(), "/Output", "/Fish_Behaviour_Data.csv"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                 Exploring                  #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# ---------------------------------------------------------------------------- #
#Function for visualizing the 'variable' for each dose with a geom_smooth, faceted by chemical... each geom_smooth represents the mean of 9 fish from that dose group. All 50 minutes of the experiment by default
# ---------------------------------------------------------------------------- #

gg_smooth_AllChems <- function (data, y, ylab, min_start = 0, min_end = 50) {
  start <- c(20,30,40)
  end <- c(25,35,45)
  dark <- data.frame(start, end)
  data %>%
    filter(time_end >= min_start & time_end <= min_end) %>%
    group_by(plate_id, Dose_mg_L, time_end) %>%
    mutate_(meanvalue = interp(~ mean(na.omit(y)), y = as.name(y))) %>% # Must use lazy evaluation in the pipe and some depracated form of mutate_
    summarise(meanvalue = unique(meanvalue)) %>%
    ggplot(data = ., mapping = aes(x = time_end, y = meanvalue, group = Dose_mg_L)) +
    geom_point(size = 0.1) +
    geom_smooth(se = FALSE, aes(color = as.factor(Dose_mg_L))) +
    scale_colour_viridis_d() +
    geom_rect(data = dark, inherit.aes = FALSE, # Adding dark rectangles to represent the dark cycles
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = 'black', alpha = 0.2) +
    geom_vline(xintercept = 20, linetype = "dashed") +
    facet_wrap(~plate_id, ncol = 6, strip.position = "top") +
    #theme_classic() +
    theme(strip.background = element_blank(), panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab("Time [min]") +
    ylab(ylab) +
    xlim(min_start, min_end) +
    guides(color = guide_legend(title = "Dose (mg/L)"))
}


#Total activity Distance
gg_smooth_AllChems(data = fishBehavDat, y = "totaldist", ylab = "Total activity distance")
gg_smooth_AllChems(data = fishBehavDat, y = "totaldist", ylab = "Total activity distance", min_start = 20)
#Interesting chemicals to explore individually... 1octanol, 24DMP, 24DNP, 4TPP, Aniline, BPA, Clofibric Acid, Dimethylformamide, DMSO, EE2, Fadrozole, Fluoxetine, Flutamide, Prochloraz, Vinclozolin... just by briefly scanning this output
#Large activity Distance
gg_smooth_AllChems(data = fishBehavDat, y = "lardist", ylab = "Large activity distance")
#Interesting chemicals to explore individually... 1octanol, 24DMP, 24DNP, 4TOP, 4TPP, BPA, Clofibric Acid, Dimethylformamide, DMSO, EE2, Fadrozole, Fluoxetine, Flutamide, Glyphosate, Prochloraz, SDS, TGSH, Vinclozolin... just by briefly scanning this output
#Small activity Distance
gg_smooth_AllChems(data = fishBehavDat, y = "smldist", ylab = "Small activity distance")
#Interesting chemicals to explore individually... 1octanol, 24DMP, 24DNP, 34DCA, 4TPP, Anailine, BPA, BPAF, Clofibric Acid, DES, Dimethylformamide, DMSO, EE2, Fadrozole, Fluoxetine, Flutamide, Prochloraz, SDS, TEG, TGSH, Vinclozolin... just by briefly scanning this output
#Active Duration
gg_smooth_AllChems(data = fishBehavDat, y = "activedur", ylab = "Total Activity Duration [0-60 seconds]")
#Interesting chemicals to explore individually... 17beta_estradiol, 24DNP, 4TOP, Clofibric Acid, Dimethylformamide, DMSO, Flutamide, Prochloraz, TEG... just by briefly scanning this output
#Large Duration
gg_smooth_AllChems(data = fishBehavDat, y = "lardur", ylab = "Total Activity Duration [0-60 seconds]")
#Small Duration
gg_smooth_AllChems(data = fishBehavDat, y = "smldur", ylab = "Total Activity Duration [0-60 seconds]")
#Inactive Duration
gg_smooth_AllChems(data = fishBehavDat, y = "inadur", ylab = "Total Activity Duration [0-60 seconds]")
# Counts... not very informative


# ---------------------------------------------------------------------------- #
#And to get a better look at what is driving the observed effect, lets get a higher resolution image... for just one chemical
# ---------------------------------------------------------------------------- #

smooth_oneChem <- function (data, chemical, y, ylab, min_start = 0, min_end = 50) {
  start <- c(20,30,40)
  end <- c(25,35,45)
  dark <- data.frame(start, end)
  data %>%
    filter(time_end >= min_start & time_end <= min_end) %>%
    filter(plate_id == chemical) %>%
    group_by(plate_id, Dose_mg_L, time_end) %>%
    mutate_(meanvalue = interp(~ mean(na.omit(y)), y = as.name(y))) %>%
    ungroup() %>%
    group_by(Dose_mg_L, replicate) %>%
    ggplot(aes(x = time_end, y = .[[y]], group = as.factor(replicate), color = as.factor(replicate))) +
    geom_line() +
    scale_colour_viridis_d() +
    geom_smooth(aes(x = time_end, y = meanvalue), se = FALSE, color = "black") +
    geom_rect(data = dark, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = 'black', alpha = 0.2) +
    facet_wrap(~Dose_mg_L, ncol = 3) +
    geom_vline(xintercept = 20, linetype = "dashed") +
    theme(strip.background = element_blank(), panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab("Time [min]") +
    ylab(ylab) +
    labs(title = chemical) +
    xlim(min_start, min_end) +
    theme(legend.position = "none")
}

# ---------------------------------------------------------------------------- #
#Interesting chemicals to explore individually for totaldist... 1octanol, 24DMP, 24DNP, 4TPP, Aniline, BPA, Clofibric Acid, Dimethylformamide, DMSO, EE2, Fadrozole, Fluoxetine, Flutamide, Prochloraz, Vinclozolin...
#Total Distance
# ---------------------------------------------------------------------------- #
smooth_oneChem(data = fishBehavDat, chemical = "1octanol", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "1octanol", y = "totaldist", ylab = "Total activity distance", min_start = 20)
smooth_oneChem(data = fishBehavDat, chemical = "24DMP", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "24DNP", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "4TPP", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "34DCA", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "Aniline", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "BPA", y = "totaldist", ylab = "Total activity distance") #Seems that one fish is driving most of the effect. Probably best to remove these outliers and replace as NA
smooth_oneChem(data = fishBehavDat, chemical = "Clofibric Acid", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "Dimethylformamide", y = "totaldist", ylab = "Total activity distance") # The highest dose seems to be driving some sort of effect for a few animals
smooth_oneChem(data = fishBehavDat, chemical = "DMSO", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "EE2", y = "totaldist", ylab = "Total activity distance") # An outlier in the third dose
smooth_oneChem(data = fishBehavDat, chemical = "Fadrozole", y = "totaldist", ylab = "Total activity distance") # A clearer pattern of effect of chemical dose on movement distance
smooth_oneChem(data = fishBehavDat, chemical = "Fluoxetine", y = "totaldist", ylab = "Total activity distance") # More effects highest dose
smooth_oneChem(data = fishBehavDat, chemical = "Flutamide", y = "totaldist", ylab = "Total activity distance") # Highest dose increasing movement distance for a few fish
smooth_oneChem(data = fishBehavDat, chemical = "Prochloraz", y = "totaldist", ylab = "Total activity distance")
smooth_oneChem(data = fishBehavDat, chemical = "Vinclozolin", y = "totaldist", ylab = "Total activity distance") # A few individuals driving these effects in the highest dose
#Conclusion, more replicates per dose group would have been beneficial for a higher throughput analysis of the effects of chemicals on the swim behaviour of fish at different doses. The next step to first convert the raw behavioural data into similarity scores, then see if there are any statistically significant effects of chemicals on the swim behaviour of fish (ANCOVA), and then a post-hoc DUnnett's test to see what doses are statistically different fromt he control dose group.



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                 Analyzing                  #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#NOTE: Paramaterize this section to analyze different endpoints (similiarity scores vs area under the curve) and different variables (total distance, active duration, inactive counts etc.)

lfishBehavDat <- fishBehavDat %>%
  mutate(value = totaldist, dose = Dose_mg_L) %>% #We will specify total activity distance
  filter(time_end > 20) %>% # For the analysis, we don't want the 20 minute acclimation period
  mutate(time_end = time_end-20) %>% # Display time 1-30 minutes
  select(plate_id, embryo_id, is_VC, time_end, value) %>%
  na.omit() # Pipeline is picky about the data frame

lfishBehavDat <- split(as.data.frame(lfishBehavDat), ~ plate_id)

#Begin Analysis pipeline adapted from Hsieh et al 2018
##Calculate Similarity Endpoints###
#This is a single endpoint for all 50 observations of a given animal...
#Create endpoints based on similarity of movement patterns in time segments
simi_endps <- list("ld_pearson" = seq(1, 30, by = 1)) #One per minute
simi_endpoints <- list()
for(i in chemicalNames){
  temp <- create_simi_endpoints(lfishBehavDat[[i]], segments = simi_endps, metric = "pearson")
  simi_endpoints[[i]] <- temp[[1]]
  rm(temp)
}
simi_endpoints

##Normalize Data##
#Normalizing the data to the internal control groups
simi_norm <- lapply(simi_endpoints, simi_normalize)
simi_norm

## Specify dose... again... Will have to play with this later##
#Specify doses, then also add a Group Variable to use for an ANCOVA... Then transform the dose column to a factor
 # for(i in chemicalNames){
 #   simi_norm[[i]] <- dose_replacement(x = simi_norm[[i]], Highdose = HighDose[i]) %>%
 #     group_by(dose)
 #   simi_norm[[i]][,"Group"] = rep(c("A","B","C"), times = 6, each = 3)
 # }
 # simi_norm

##Statistics##
for(i in chemicalNames){
  simi_norm[[i]]$dose <- as.factor(simi_norm[[i]]$dose)
  simi_norm[[i]][,"Group"] = factor(rep(c("A","B","C"), times = (nrow(simi_norm[[i]])/9), each = 3))
} # Adding Group column and factoring dose variable
as_tibble(simi_norm[["BPA"]])

#Obtaining Summary Statistics
summarystats_list <- lapply(simi_norm, summarystats)
summarystats_list[["BPA"]]

#ANCOVA
ancova_list <- sapply(simi_norm, combinedancova)
print(ancova_list[["BPA"]])
print(ancova_list[["4TPP"]])
#4TPP similarity scores statistically significant from each other
ancova_comb <- as_tibble(ldply(ancova_list))
ancova_comb <- ancova_comb %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = `Pr(>F)` < 0.05,
    true = TRUE,
    false = FALSE
  ))
ancova_comb
#Any significant results?
ancova_comb %>%
  filter(is.significant == TRUE)
# Total distance - Yep, 4TPP is significant... be cautious look at the raw data
# Active Duration - Yep, DMSO is significant... be cautious look at the raw data
# Times gone Inactive (inact) - Yep, BPAF, Fenitrothion, TGSH... and fadrozole 'group' explain results


#ANOVA
anova_list <- sapply(simi_norm, combinedanova)
print(anova_list[["BPA"]])
print(anova_list[["4TPP"]])
#4TPP similarity scores statistically significant from each other
anova_comb <- as_tibble(ldply(anova_list))
anova_comb <- anova_comb %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = `Pr(>F)` < 0.05,
    true = TRUE,
    false = FALSE
  ))
anova_comb
#Any significant results?
anova_comb %>%
  filter(is.significant == TRUE)
# Total Distance - Yep, 4TPP is significant... be cautious look at the raw data
# Active Duration - Yep, DMSO is significant... be cautious look at the raw data
# Times gone Inactive (inact) - Yep, BPAF, Fenitrothion is significant... be cautious look at the raw data


#Post Hoc Test - Dunnett's test
# combineddunnett <- function(df){
#   df$dose = factor(str_split(string = df$embryo_id, pattern = "_", simplify = TRUE)[,1])
#   DunnettTest(x = df$endpoint_value, g = df$dose)
# }
dunnett_list <- sapply(simi_norm, combineddunnett) %>%
  setNames(., chemicalNames) %>% 
  as.array()
print(dunnett_list[["4TPP"]])
#Nothing
#Dunnett's Test requires a bunch of wrangling
dunnett_temp <- list() #Take Dunnett's test results without any of the fancy summary information and shove it into a named list
for (i in chemicalNames) {
  dunnett_temp[[names(dunnett_list[i])]] <-
    dunnett_list[[i]] %>% as.data.frame() %>% rownames_to_column(var = "dose") #Coerce to a data frame temporarily so what we can take the row names of the reults and turn them into a variable with rownames_to_column
}
dunnett_comb <- ldply(dunnett_temp) #this function combines all of the lists together and gives them a variable name according to the chemical
dunnett_comb$dose = substr(dunnett_comb$dose,
                           start = 1,
                           stop = nchar(dunnett_comb$dose) - 6) #Here we are fixing the dose column... the dose column has the test dose related to the control... but we just want to see what the test dose is without it giving us redundant information about the comparison to the control for every observation...
dunnett_comb
dunnett_comb <- as_tibble(dunnett_comb)
dunnett_comb <- dunnett_comb %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = pval < 0.05,
    true = TRUE,
    false = FALSE
  ))
dunnett_comb
# Any significant results?
dunnett_comb %>%
  filter(is.significant == TRUE)
# Total Distance - Nope
# Active Duration - DMSO - Dose 4 and Dose 5
# Times fish went inactive during experiment - BPAF - Dose 2, TGSH - Dose 4


###Export Data###
## Total Distance Value/Variable
for(i in chemicalNames){
  write.csv(simi_norm[[i]], file = paste0("Output/", i,"/",i,"_simi_norm.csv"))
}
for(i in chemicalNames){
  write.csv(anova_list[[i]], file = paste0("Output/", i,"/",i,"_anova.csv"))
}
for(i in chemicalNames){
  write.csv(dunnett_list[[i]], file = paste0("Output/", i,"/",i,"_dunnett.csv"))
}
for(i in chemicalNames){
  write.csv(ancova_list[[i]], file = paste0("Output/", i,"/",i,"_ancova.csv"))
}

write.csv(x = dunnett_comb, file = "Output/Dunnett_Results_Simi_endpoint.csv")
write.csv(x = ancova_comb, file = "Output/ANCOVA_Results_Simi_endpoint.csv")
write.csv(x = anova_comb, file = "Output/ANOVA_Results_Simi_endpoint.csv")

#Export Normalized Similarity scores in a tibble
doseData <- doseData %>%
  mutate(dose = Dose) #Renaming dose variable to inner_join
simi_norm_tib <- as_tibble(rbindlist(simi_norm)) #Convert from list to tibble
simi_norm_tib <- simi_norm_tib %>%
  inner_join(doseData) %>%
  mutate(is_VC = as.integer(is_VC)) %>%
  select(plate_id, dose, Dose_mg_L, is_VC, Group, endpoint, endpoint_value_norm)
simi_norm_tib

write_csv(simi_norm_tib, file = "Output/simi_norm_data.csv")

## Active Duration Value/Variable

## Amount of times fish went inactive during the experiment (inact) Value/Variable

####Plotting####
# behaviourplots <- list()
# for(i in chemicalNames) {
#   behaviourplots[[i]] <-
#     simi_norm[[i]] %>%
#     ggplot(aes(x = as.character(dose), y = endpoint_value_norm, group = dose)) +
#     geom_boxplot(outlier.shape = NA, width = 0.5) +
#     geom_jitter(position = position_jitter(width = 0.2,
#                                            height = 0, seed = 42069),
#                 colour = "black") +
#     labs(title = paste0(i)) + # x = "Dose (µg/L)", y = "Response") +
#     theme_classic() +
#     theme(axis.title.x = element_blank(), axis.title.y = element_blank())
# }
# print(behaviourplots[[8]])
# 
# multiplot(plotlist = behaviourplots, cols = 5 )
# 
# #Playing with plots
# behaviourplots <- list()
# for(i in chemnames) {
#   behaviourplots[[i]] <-
#     simi_norm[[i]] %>%
#     ggplot(aes(x = as.numeric(dose), y = endpoint_value_norm)) +
#     geom_boxplot(aes(group = dose), outlier.shape = NA, width = 0.5) +
#     geom_smooth(se = FALSE, method = "auto") +
#     geom_jitter(width = 0.2,
#                 height = 0,
#                 colour = "black") +
#     labs(title = paste0(i), x = "Dose (µg/L)", y = "Response") +
#     scale_x_continuous(breaks = c(1:6), labels = levels(simi_norm[[i]]$dose)) +
#     theme_classic()
# }
# print(behaviourplots[[1]])

#Plotting from tibble
simi_norm_tib %>%
  filter(plate_id == "DMSO") %>%
  group_by(plate_id, Dose_mg_L, Group) %>%
  summarise(Dose_mg_L = Dose_mg_L, endpoint_value_norm = endpoint_value_norm, Group = Group) %>%
  ggplot(aes(x = as.factor(Dose_mg_L), y = endpoint_value_norm, group = Dose_mg_L)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
      # geom_jitter(position = position_jitter(width = 0.2,
      #                                        height = 0, seed = 42069),
      #             colour = "black") +
  geom_jitter((aes(colour = Group)), position = position_jitter(width = 0.2, height = 0, seed = 42069)) +
  labs(x = "Dose (mg/L)", y = "Response") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~plate_id)



############ NOT DONE ###############
####Setup Data for RCurvep####
simi_norm_tib_4rcurvep <- simi_norm_tib %>%
  rename(resp = endpoint_value_norm,
         chemical = plate_id) %>%
  filter(is_VC == 0) %>%
  mutate(conc = log10(Dose_mg_L)) %>%
  select(endpoint, chemical, conc, resp)
glimpse(simi_norm_tib_4rcurvep)

#Estimate BMR - bootstrap - simulate curves
set.seed(42069)

simi_norm_tib_4rcurvep_act <- combi_run_rcurvep(
  simi_norm_tib_4rcurvep, 
  n_samples = 100, 
  keep_sets = c("act_set"), 
  TRSH = seq(5, 95, by = 5)
)

simi_norm_tib_4rcurvep_act[["result"]][["act_set"]] %>%
  ggplot(aes(x = wConc, y = wResp, group = sample_id)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~chemical, scales = "free")

simi_norm_tib_4rcurvep_act_Summary_Output <- summarize_rcurvep_output(simi_norm_tib_4rcurvep_act)

bmr_output <- estimate_dataset_bmr(simi_norm_tib_4rcurvep_act, plot = TRUE)
bmr_output

plot(bmr_output)

#Calculating the BMD based off of the BMR using the 'Hill' model
fit_dat <- run_fit(create_dataset(simi_norm_tib_4rcurvep), hill_pdir = 1, n_samples = 10, modls = "hill")
fit_dat

fit_dat_summary_output <-
  summarize_fit_output(fit_dat, thr_resp = bmr_output$outcome$bmr_ori)
fit_dat_summary_output

fit_dat_summary_output[["act_summary"]] %>% View()

#Calculating the BMD based off of the BMR using the 'CurveP' model





##Testing running for each individual chemical
test <- split(x = simi_norm_tib_4rcurvep, f = ~chemical, drop = FALSE) 


#Estimating the BMR - Bootstrap
test_act <- combi_run_rcurvep(
  test[[1]], 
  n_samples = 100, 
  keep_sets = c("act_set"), 
  TRSH = seq(5, 95, by = 5), # test all candidates, 5 to 95
  RNGE = 10
)
#summarize_rcurvep_output(test_act)
test_act[["result"]][["act_set"]] %>%
  ggplot(aes(x = wConc, y = wResp, group = sample_id)) + 
  geom_point() +
  geom_line()

bmr_output <- estimate_dataset_bmr(test_act, plot = TRUE)
bmr_output$outcome

#Calculating the BMD based off the BMR
fit_dat <- run_fit(create_dataset(simi_norm_tib_4rcurvep), hill_pdir = 1, n_samples = 10, modls = "hill")
fit_dat

fit_dat_summary_output <-
  summarize_fit_output(fit_dat, thr_resp = bmr_output$outcome$bmr_ori)
fit_dat_summary_output

fit_dat_summary_output[["act_summary"]] %>% filter(hit_confidence >= 0.1) %>% View()