#### Creating your data frame ####
#Load all the required packaged
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

#Importing the data
data <- fread("~/MolToxBackup/Behaviour/BPA/BPA_(copy).XLS",
              sep = "auto",
              data.table = FALSE, #creates data.frame
              drop = c(1, 3, 8), #drop unnecessary columns
              blank.lines.skip = TRUE 
)
#More Cleaning
data <- data %>%
  drop_na() %>% #delete all empty wells
  filter(row_number() %% 2 != 0) #filter and remove every even row. For x != y, x= how often to delete and y= where to begin

#Create total/summary columns
data <- data %>%
  mutate(tmovect = smlct + larct) %>%
  mutate(tmovedur = smldur + lardur) %>%
  mutate(tmovedist = smldist + lardist)

#### Binning by minute (60s) ####
breaks <- seq(from = 0, to = 3060, by = 60)
tags <- c("[0-60)", "[60-120)", "[120-180)", "[180-240)","[240-300)", "[300-360)", "[360-420)",
          "[420-480)", "[480-540)", "[540-600)", "[600-660)", "[660-720)", "[720-780)",
          "[780-840)", "[840-900)", "[900-960)", "[960-1020)", "[1020-1080)",
          "[1080-1140)", "[1140-1200)", "[1200-1260])","[1260-1320)", "[1320-1380)", "[1380-1440)", "[1440-1500)",
          "[1500-1560)", "[1560-1620)", "[1620-1680)", "[1680-1740)","[1740-1800)", "[1800-1860)","[1860-1920)", "[1920-1980)",
          "[1980-2040)", "[2040-2100)", "[2100-2160)", "[2160-2220)", "[2220-2280)",
          "[2280-2340)", "[2340-2400)", "[2400-2460)", "[2460-2520)", "[2520-2580)",
          "[2580-2640)", "[2640-2700)", "[2700-2760)", "[2760-2820)", "[2820-2880)",
          "[2880-2940)", "[2940-3000)", "[3000-3060)")


# Creating objects for each dose group
#### Dose 1 ####
Dose1 <- data[data$Treatment == str_extract(data$Treatment, "Dose1.") ,]
Dose1 <- Dose1[!is.na(Dose1$Treatment),]
#Binning
group_tags1 <- cut(Dose1$start,
                  breaks = breaks,
                  include.lowest = TRUE,
                  right = FALSE,
                  labels = tags)
minute_groups1 <- factor(group_tags1, levels = tags, ordered = TRUE)
#Bin w Mean replicates
MDose1 <- Dose1 %>%
  aggregate(by = list(group_tags1), mean) %>%
  subset(select = -Treatment)

#### Dose 2 ####
Dose2 <- data[data$Treatment == str_extract(data$Treatment, "Dose2.") ,]
Dose2 <- Dose2[!is.na(Dose2$Treatment),]
#Binning
group_tags2 <- cut(Dose2$start,
                  breaks = breaks,
                  include.lowest = TRUE,
                  right = FALSE,
                  labels = tags)
minute_groups2 <- factor(group_tags2, levels = tags, ordered = TRUE)
#Bin w Mean replicates
MDose2 <- Dose2 %>%
  aggregate(by = list(group_tags2), mean) %>%
  subset(select = -Treatment)

#### Dose 3 ####
Dose3 <- data[data$Treatment == str_extract(data$Treatment, "Dose3.") ,]
Dose3 <- Dose3[!is.na(Dose3$Treatment),]

#Binning
group_tags3 <- cut(Dose3$start,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = tags)
minute_groups3 <- factor(group_tags3, levels = tags, ordered = TRUE)
#Bin w Mean replicates
MDose3 <- Dose3 %>%
  aggregate(by = list(group_tags3), mean) %>%
  subset(select = -Treatment)

#### Dose 4 ####
Dose4 <- data[data$Treatment == str_extract(data$Treatment, "Dose4.") ,]
Dose4 <- Dose4[!is.na(Dose4$Treatment),]

#Binning
group_tags4 <- cut(Dose4$start,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = tags)
minute_groups4 <- factor(group_tags4, levels = tags, ordered = TRUE)
#Bin w Mean of replicates and create new object
MDose4 <- Dose4 %>%
  aggregate(by = list(group_tags4), mean) %>%
  subset(select = -Treatment)

#### Dose 5 ####
Dose5 <- data[data$Treatment == str_extract(data$Treatment, "Dose5.") ,]
Dose5 <- Dose5[!is.na(Dose5$Treatment),]

#Binning
group_tags5 <- cut(Dose5$start,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = tags)
minute_groups5 <- factor(group_tags5, levels = tags, ordered = TRUE)
#Bin w Mean of replicates and create new object
MDose5 <- Dose5 %>%
  aggregate(by = list(group_tags5), mean) %>%
  subset(select = -Treatment)

#### Dose 6 ####
CC <- data[data$Treatment == str_extract(data$Treatment, "Control.") ,]
CC <- CC[!is.na(CC$Treatment),]

#Binning
group_tagsc <- cut(CC$start,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = tags)
minute_groupsc <- factor(group_tagsc, levels = tags, ordered = TRUE)
#Bin w Mean of replicates and create new object
MCC <- CC %>%
  aggregate(by = list(group_tagsc), mean) %>%
  subset(select = -Treatment)

#Plotting
colors <- c("Dose 1" = "navy","Dose 2" = "dodgerblue","Dose 3" = "lightskyblue","Dose 4" = "lightsteelblue1","Dose 5" = "snow","CC" = "black")
#### Total Movement Count ####
Mcount = ggplot() +
  geom_line(data = MDose1, aes(x = start/60, y = tmovect, color = "Dose 1")) + 
  geom_line(data = MDose2, aes(x = start/60, y = tmovect, color = "Dose 2")) +
  geom_line(data = MDose3, aes(x = start/60, y = tmovect, color = "Dose 3")) +
  geom_line(data = MDose4, aes(x = start/60, y = tmovect, color = "Dose 4")) +
  geom_line(data = MDose5, aes(x = start/60, y = tmovect, color = "Dose 5")) +
  geom_line(data = MCC, aes(x = start/60, y = tmovect, color = "CC")) +
  labs(x = "Time (m)", y = "Total Movement Counts", color = "Legend") +
  ggtitle("Total larval movement counts") +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = 20, linetype = "dotted")
print(Mcount)

#### Total Movement Duration ####
Mdur = ggplot() +
  geom_line(data = MDose1, aes(x = start/60, y = tmovedur, color = "Dose 1")) + 
  geom_line(data = MDose2, aes(x = start/60, y = tmovedur, color = "Dose 2")) +
  geom_line(data = MDose3, aes(x = start/60, y = tmovedur, color = "Dose 3")) +
  geom_line(data = MDose4, aes(x = start/60, y = tmovedur, color = "Dose 4")) +
  geom_line(data = MDose5, aes(x = start/60, y = tmovedur, color = "Dose 5")) +
  geom_line(data = MCC, aes(x = start/60, y = tmovedur, color = "CC")) +
  labs(x = "Time (m)", y = "Total Movement duration", color = "Legend") +
  ggtitle("Duration of larval movements") +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = 20, linetype = "dotted")
print(Mdur)

#### Total Movement Distance ####
Mdist = ggplot() +
  geom_line(data = MDose1, aes(x = start/60, y = tmovedur, color = "Dose 1")) + 
  geom_line(data = MDose2, aes(x = start/60, y = tmovedur, color = "Dose 2")) +
  geom_line(data = MDose3, aes(x = start/60, y = tmovedur, color = "Dose 3")) +
  geom_line(data = MDose4, aes(x = start/60, y = tmovedur, color = "Dose 4")) +
  geom_line(data = MDose5, aes(x = start/60, y = tmovedur, color = "Dose 5")) +
  geom_line(data = MCC, aes(x = start/60, y = tmovedur, color = "CC")) +
  labs(x = "Time (m)", y = "Total Movement distance", color = "Legend") +
  ggtitle("Distance of larval movements") +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = 20, linetype = "dotted")
print(Mdist)

#### Totals ####
