####Libraries####
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)

####Directory####

#To avoid errors on different devices...
#Convert .XLS files to .csv files
cd ~/MolToxLab/Behaviour/Data/All #change directory to the 'All' folder
ls #list your files, take a look these should be your files
for i in *.XLS
do soffice --headless --convert-to csv "$i"
done #Run a loop to convert all the .XLS files to csv
#I will post these .csv files onto the git repository so you don't have to do this

getwd() #Should be something like ~/MolToxLab/Behaviour if you open script from file
folderNames <- list.files("Data/Sub")
folder <- folderNames
filesinfolder <- list.files(path = paste0("Data/Sub/", folder, "/"), pattern = "*.csv")

####Importing and Cleaning####

#Import as list
list <- list()
for (i in 1:length(filesinfolder)) {
  list[[i]] <- read_delim(file = paste0("Data/All/", filesinfolder[[i]]),
                          delim = "\t",
                          col_names = TRUE,
                          col_types = "cclnninninninnin",
                          col_select = -c("an" ,"inadist", "emptyct", "emptydur"),
                          na = c("", 'NA', "NA", "\t NA", "\tNA")
  )
}
names(list) <- folder #Name List objects
#Delete every duplicate row
list2 <- lapply(list, FUN = distinct)
#Get rid of Empty wells ("NA" values)
list2 <- lapply(list2, FUN = drop_na)
#Drop all data past 50 minutes
only50mins <- function(x) {
  x <- x[-c(2701:2754), ]
}
list2 <- lapply(list2, FUN = only50mins)
#Create one big data frame... now you can use it in tidyverse :D
data <- rbindlist(list2, idcol = "Chemical")
rm(list)
#Spread out 'Treatment' Column into 'Dose' and 'Group'
data2 <- data %>%
  mutate(Treatment = rep(c("Dose1_A", "Dose1_B", "Dose1_C","Dose2_A", "Dose2_B", 
                         "Dose2_C", "Dose3_A", "Dose3_B", "Dose3_C", "Dose4_A", 
                         "Dose4_B", "Dose4_C", "Dose5_A", "Dose5_B", "Dose5_C", 
                         "Dose6_A", "Dose6_B", "Dose6_C"), 
                         times = 1300, each = 3)
  ) %>%
  separate(col = Treatment, into = c('Dose', 'Group'), convert = TRUE)

#Import Dose list
AllDoses <- read.csv(file  = "Highest_Dose_Aug_1.csv", skip = 1, header = TRUE)
AllDoses2 <- AllDoses %>%
  gather(key = Dose, value = "Dose(mg/L)", Dose1:Dose6)

#Create the final data frame
FinalData <- data2 %>%
  inner_join(AllDoses2) %>%
  select(-Dose)

rm(AllDoses, data, data2)

####Analysis####

#####Plots and Visuals####
start <- c(20,30,40)
end <- c(25,35,45)
dark <- data.frame(start, end)
#Inactivity Counts
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avginact = mean(inact)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avginact) %>%
  ggplot(aes(x = time, y = avginact)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Inactivity count")

#Inactivity Duration
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avginadur = mean(inadur)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avginadur) %>%
  ggplot(aes(x = time, y = avginadur)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Inactivity Duration")

#Small Activity Counts
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgsmlct = mean(smlct)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgsmlct) %>%
  ggplot(aes(x = time, y = avgsmlct)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Small activity count")

#Small Activity Duration
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgsmldur = mean(smldur)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgsmldur) %>%
  ggplot(aes(x = time, y = avgsmldur)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Small activity duration")

#Small Activity Distance
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgsmldist = mean(smldist)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgsmldist) %>%
  ggplot(aes(x = time, y = avgsmldist)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Small activity distance")

#Large Activity Counts
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avglarct = mean(larct)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avglarct) %>%
  ggplot(aes(x = time, y = avglarct)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Large activity count")

#Large Activity Duration
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avglardur = mean(lardur)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avglardur) %>%
  ggplot(aes(x = time, y = avglardur)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Large activity duration")

#Large Activity Distance
FinalData %>%
  mutate(time = end/60) %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avglardist = mean(lardist)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avglardist) %>%
  ggplot(aes(x = time, y = avglardist)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Large activity distance")
