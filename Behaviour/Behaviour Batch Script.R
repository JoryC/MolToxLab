#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####                 Version Control            #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#R version 4.1.1 (2021-08-10) -- "Kick things"
#Platform: x86_64-pc-linux-gnu (64-bit) 5.11.0-34-generic / Ubuntu 20.04.03
#Desktop: GNOME 3.36.5
#Hardware: CPU - Intel Core i5-9400F 6 core 4.1GHz / RAM - 15924MiB
#ggplot2 Package: Version 3.3.5
#dplyr Package: Version 1.0.7
#tidyr: Version 1.1.3
#readr Package: Version 2.0.0
#stringr: 1.4.0
#data.table: 1.14.0
#car: 3.0-11

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                   Libraries                ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(car)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####                Directory                    ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####            Importing and Cleaning          #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#Scroll to the bottom of this section to easily import the data from FinalData.csv
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
  select(-Dose) %>%
  mutate(totaldist = smldist + lardist,
         totaldur = smldur + lardur,
         totalct = smlct + larct,
         time = end/60,
         velocity = totaldist/time)

rm(AllDoses, data, data2)

#Write the finished data frame as a .csv for others to easily import
write.csv(FinalData, "~/MolToxLab/Behaviour/FinalData.csv")

#The lazy way to import all the data
FinalData <- read_csv("FinalData.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####                Analysis                 #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Comments and unused code ####
#Bin 1 - 5 minutes before first dark cycle
#baselinebin <- FinalData %>%
#  filter(time == 15:20)
#Bin 2 - Three 5 minute dark cycle chunks - total 15 mins dark
#darkbin <- FinalData %>%
#  filter(time == 20:25 & 30:35 & 40:45)
#Bin 3 - Three 5 minute light cycle chunks - total 15 mins light
#lightbin <- FinalData %>%
#  filter(time == 25:30 & 35:40 & 45:50)
#All bins
#bins <- as.data.frame(c(baselinebin, darkbin, lightbin))
#### Pre-Analysis ####
#Let's run an two-way ANOVA on total activity distance
#Independent variables: Dose and time
#Dependent: In this case, total activity distance

#### All Chemicals totaldistance Boxplot ####
#Let's boxplot out Dose on the x axis, totaldistance on the y axis, and facet wrap by chemical...
FinalData %>%
  group_by(Chemical, `Dose(mg/L)`) %>%
  select(Chemical, `Dose(mg/L)`, totaldist) %>%
  ggplot() + 
  geom_boxplot(aes(group = `Dose(mg/L)`, x = as.factor(`Dose(mg/L)`), y = totaldist)) +
  facet_wrap(~Chemical, scales  = "free")

#### 17beta_estradiol total distance summary ####
#And now let's take a deeper dive into just one chemical...
#First filter by chemical (17beta_estradiol in this case),
#Check out the means and sd's

#### 17beta_estradiol total distance boxplot ####
#Defining our light and dark time periods for the geom_rect plot coming up
start <- c(20,30,40)
end <- c(25,35,45)
dark <- data.frame(start, end)
#Creating our ggplot pipeline per chemical (we are using 17bestradiol here)
#We're plotting out 'totaldistance' over 'time' and facet wrapping by the 'dose'
tesboxplot <- FinalData %>%
  filter(Chemical == "17beta_estradiol", time >= 5) %>%
  group_by(time, `Dose(mg/L)`) %>%
  select(totaldist, time, `Dose(mg/L)`) %>%
  arrange(`Dose(mg/L)`) %>%
  ggplot() + 
  geom_boxplot(aes(group = time, x = time, 
                   y = totaldist)) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~`Dose(mg/L)`, scales = "free") +
  xlim(10, 50)
  
ggsave(filename = "testboxplot.pdf", plot = tesboxplot, 
       path = paste0(getwd(), "/Images/Boxplots"), 
       device = "pdf",
       width = 1920,
       height = 1080,
       units = "px",
       scale = 3)

#### 17beta_estradiol ANOVA ####
#Now that we've done a bit of visualizing, let's try to model this
#First let's define the data we want to model... here we'll try to model 17betaestradiol
b17estradiol <- FinalData %>%
  filter(Chemical == "17beta_estradiol", time >= 5) %>%
  group_by(time, `Dose(mg/L)`) %>%
  select(totaldist, time, `Dose(mg/L)`) %>%
  arrange(`Dose(mg/L)`)
#Now we are ready to plug the data into the ANOVA model
totaldistmodel <- aov(totaldist ~ as.factor(`Dose(mg/L)`) * as.factor(time), data = b17estradiol)
summary(totaldistmodel) #View the results
#Model residuals subset
resid <- totaldistmodel$residuals
#Histogram of residuals (normality?) Should look like a normal dist
hist(resid, xlim = range(-1500:1500)) #Is it normally distributed?

#### Assumptions: 17beta_estradiol Levene's test equal variance & Normality ####
#Let's test that the variances are approximately equal
test <- leveneTest(totaldist ~ as.factor(`Dose(mg/L)`) * as.factor(time), data = b17estradiol)
#If below p0.05 then it's not equal variance
if (test$`Pr(>F)`[1] >= 0.05) {
  print("Equal Variance")
} else {
  print("Not Equal Variance")
}

#### 17beta_estradiol Tukey's test ####
#Great, so now let's try and see which treatment groups differ from one another
#To do this we will do a Tukey test for multiple comparisons
tukey <- TukeyHSD(totaldistmodel, conf.level=.99)
#Which dose groups are different from one another?
which(tukey$'as.factor(`Dose(mg/L)`)'[,"p adj"] <= 0.01)
#Which times binned by minute are different from each other?
#Note: we expect high significance between different bins (light-dark, light-baseline, dark-baseline, dark-light, etc.)
#Note 2: Is this really that informative at all?
which(tukey$`as.factor(time)`[,"p adj"] <= 0.01)
#Which were significant taking into account both time and dose
#Note: this will be a long list
#How to summarise interpret this data? 
which(tukey$'as.factor(`Dose(mg/L)`):as.factor(time)'[,"p adj"] <= 0.01)



#### Analysis - All Chemicals ####
#And to do this of all chemicals we can create a loop to do whatever function we want
#### Summaries ####
#
SummaryData <- FinalData %>%
    filter(time >= 5) %>%
    group_by(Chemical, `Dose(mg/L)`, time) %>%
    summarise(meaninact = mean(inact), sdinact = sd(inact),
              meaninadur = mean(inadur), sdinadur = sd(inadur),
              meansmlct = mean(smlct), sdsmlct = sd(smlct),
              meansmldur = mean(smldur), sdsmldur = sd(smldur),
              meansmldist = mean(smldist), sdsmldist = sd(smldist),
              meanlarct = mean(larct), sdlarct = sd(larct),
              meanlardur = mean(lardur), sdlardur = sd(lardur),
              meanlardist = mean(lardist), sdlardist = sd(lardist),
              meantotaldist = mean(totaldist), sdtotaldist = sd(totaldist),
              meantotaldur = mean(totaldur), sdtotaldur = sd(totaldur),
              meantotalct = mean(totalct), sdtotalct = sd(totalct),
              meanvelocity = mean(velocity), sdvelocity = sd(velocity)
              )
#### Boxplot Visualizations ####
#Defining variable that will be used in the function
start <- c(20,30,40)
end <- c(25,35,45)
dark <- data.frame(start, end)
#Creating our list of first-pass boxplot visualizations facet wrapped by dose for each chemical
#List of all the treatment strings
treatments <- names(FinalData[-c(1, 2, 3, 4, 5, 6, 15, 19)])
#Function to plot by treatment string
#Note: Must find a conclusion on how to define outliers
customggboxplot <- function(y, chem) {
  FinalData %>%
    filter(Chemical == chem, time >= 5) %>%
    group_by(time, `Dose(mg/L)`) %>%
    select(y, time, `Dose(mg/L)`) %>%
    arrange(`Dose(mg/L)`) %>%
    ggplot() + 
    geom_boxplot(aes_string(group = "time", x = "time", 
                     y = y)) +
    geom_rect(data = dark, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = 'black', alpha = 0.2) +
    facet_wrap(~`Dose(mg/L)`, scales = "free") +
    xlim(4.5, 50.5)
}
#Set of all the custom ggplot images
boxplotlist <- list()
for (j in 1:length(folderNames)) {
boxplotlist[[folderNames[j]]] <- sapply(X = treatments, FUN = customggboxplot, chem = folderNames[j], simplify = FALSE)
}
#Great now that we can visualize our data very simply just by typing out the paramaters we want in the format of...
#boxplotlist$chemical$treatment... i.e.
boxplotlist$BPA$lardist

#Saving Each plot locally 
#This takes a long time, don't bother running it if u dont want ~312 pdf files stored locally on your computer
for (k in 1:length(boxplotlist)) {
  for (l in 1:length(boxplotlist[[k]])) {
    ggsave(filename = paste0(names(boxplotlist[k]), "_", names(boxplotlist[[k]][l]), ".pdf"), 
           plot = boxplotlist[[k]][[l]], 
           path = paste0(getwd(), "/Images/Boxplots"), 
           device = "pdf",
           width = 1920,
           height = 1080,
           units = "px",
           scale = 3)
  }
}

#for (k in 1:length(boxplotlist[])){
#  for (l in 1:length(boxplotlist[][[l]])) {
#    ggsave(filename = paste0(names(boxplotlist[k]), "_", names(boxplotlist[][[l]][l]), ".pdf"), plot = tesboxplot, 
#           path = paste0(getwd(), "/Images/Boxplots"), 
#           device = "pdf",
#           width = 1920,
#           height = 1080,
#           units = "px",
#           scale = 3)
#  }
#}
#### ANOVAs ####
#### Normailty ####
aovdatalist <- NULL
tempaov <- NULL
for (n in 1:length(boxplotlist)){
  for (o in 1:length(boxplotlist[[n]])) {
    chem = unique(FinalData$Chemical)[n]
    var = treatments[o]
    data = FinalData %>%
      group_by(Chemical, time, `Dose(mg/L)`) %>%
      select(Chemical, var, time, `Dose(mg/L)`) %>%
      filter(Chemical == chem, time >= 5) %>%
        arrange(`Dose(mg/L)`) %>%
      ungroup()
    aovdatalist[[folderNames[n]]][[treatments[o]]] <- data
     tempaov <- aov(
        formula = as.numeric(as.matrix(aovdatalist[[n]][[o]][,2])) ~ as.factor(`Dose(mg/L)`) * as.factor(time),
        data = aovdatalist[[n]][[o]])
     resid = tempaov$residuals
     png(filename = paste0(getwd(), "/Images/Residuals Histograms/", names(aovdatalist[n]), names(aovdatalist[[n]][[o]][,2]), "residuals.png"))
     hist(resid)
     dev.off()
  }
}
#View the results
#Histogram of residuals (normality?) Should look like a normal dist
#Is it normally distributed?
aovdatalist <- NULL
tempaov <- NULL
for (n in 1:length(boxplotlist)){
  for (o in 1:length(boxplotlist[[n]])) {
    chem = unique(FinalData$Chemical)[n]
    var = treatments[o]
    data = FinalData %>%
      group_by(Chemical, time, `Dose(mg/L)`) %>%
      select(Chemical, var, time, `Dose(mg/L)`) %>%
      filter(Chemical == chem, time >= 5) %>%
      arrange(`Dose(mg/L)`) %>%
      ungroup()
    aovdatalist[[folderNames[n]]][[treatments[o]]] <- data
    tempaov <- aov(
      formula = as.numeric(as.matrix(aovdatalist[[n]][[o]][,2])) ~ as.factor(`Dose(mg/L)`) * as.factor(time),
      data = aovdatalist[[n]][[o]])
    resid = tempaov$residuals
    print(paste(names(aovdatalist[n]), names(aovdatalist[[n]][[o]][,2])))
    print(shapiro.test(x = resid))
  }
}
    
#### Equal Variance ####
#Let's test that the variances are approximately equal
aovdatalist <- NULL
tempaov <- NULL
levtest <- NULL
for (n in 1:length(boxplotlist)){
  for (o in 1:length(boxplotlist[[n]])) {
    chem = unique(FinalData$Chemical)[n]
    var = treatments[o]
    data = FinalData %>%
      group_by(Chemical, time, `Dose(mg/L)`) %>%
      select(Chemical, var, time, `Dose(mg/L)`) %>%
      filter(Chemical == chem, time >= 5) %>%
      arrange(`Dose(mg/L)`) %>%
      ungroup()
    aovdatalist[[folderNames[n]]][[treatments[o]]] <- data
    tempaov <- aov(
      formula = as.numeric(as.matrix(aovdatalist[[n]][[o]][,2])) ~ as.factor(`Dose(mg/L)`) * as.factor(time),
      data = aovdatalist[[n]][[o]])
    levtest[[folderNames[n]]][[treatments[o]]] <- leveneTest(
      y = as.numeric(as.matrix(aovdatalist[[n]][[o]][,2])) ~ as.factor(as.matrix(aovdatalist[[n]][[o]][,4])) * as.factor(as.matrix(aovdatalist[[n]][[o]][,3])),
                         data = data)
  }
}
#test <- leveneTest(totaldist ~ as.factor(`Dose(mg/L)`) * as.factor(time), data = b17estradiol)
#If below p0.05 then it's not equal variance
#leveneTest(response variable ~ group variable * , data = data)
for (n in 1:length(boxplotlist)){
  for (o in 1:length(boxplotlist[[n]])) {
      if (levtest[[n]][[o]][["Pr(>F)"]][1] >= 0.05) {
  print(paste(names(aovdatalist[n]),names(aovdatalist[[n]][[o]][,2]), "Equal Variance"))
}       else {
  next    
  #print(paste(names(aovdatalist[n]),names(aovdatalist[[n]][[o]][,2]), "Not Equal Variance"))
}
  }
}

#### Tukey's test ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####              Visualization                ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
start <- c(20,30,40)
end <- c(25,35,45)
dark <- data.frame(start, end)
#### Inactivity Counts ####
inacntsplot <- FinalData %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avginact = mean(inact)) %>%              #Normalization (mean of 3 fish per minute)
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
  ylab("Inactivity count") +
  

ggsave(filename = "Inactivity_counts.pdf", plot = inacntsplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3)

#### Inactivity Duration ####
inadurplot <- FinalData %>%
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
  ylab("Inactivity Duration") +
  xlim(5, 50)

ggsave(filename = "Inactivity_duration.pdf", plot = inadurplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Small Activity Counts ####
smactcntsplot <- FinalData %>%
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
  ylab("Small activity count") +
  xlim(5, 50)

ggsave(filename = "Small_activity_counts.pdf", plot = smactcntsplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Small Activity Duration ####
smactdurplot <- FinalData %>%
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
  ylab("Small activity duration") +
  xlim(5, 50)

ggsave(filename = "Small_activity_duration.pdf", plot = smactdurplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Small Activity Distance ####
smactdistplot <- FinalData %>%
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
  ylab("Small activity distance") +
  xlim(5, 50)

ggsave(filename = "Small_activity_distance.pdf", plot = smactdistplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Large Activity Counts ####
laractcntsplot <- FinalData %>%
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
  ylab("Large activity count") +
  xlim(5, 50)

ggsave(filename = "Large_activity_counts.pdf", plot = laractcntsplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Large Activity Duration ####
laractdurplot <- FinalData %>%
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
  ylab("Large activity duration") +
  xlim(5, 50)

ggsave(filename = "Large_activity_duration.pdf", plot = laractdurplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Large Activity Distance ####
laractdistplot <- FinalData %>%
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
  ylab("Large activity distance") +
  xlim(5, 50)

ggsave(filename = "Large_activity_distance.pdf", plot = laractdistplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Distance/Time (Velocity - sumdistance/minute) ####
velocityplot <- FinalData %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgvelocity = mean(velocity)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgvelocity) %>%
  ggplot(aes(x = time, y = avgvelocity)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  xlim(10, 50) +
  ylim(0, 150) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Velocity [dist/min]")

ggsave(filename = "Velocity.pdf", plot = velocityplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)
#### Total Activity Counts ####
totalctplot <- FinalData %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgtotalct = mean(totalct)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgtotalct) %>%
  ggplot(aes(x = time, y = avgtotalct)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Total activity count") + 
  xlim(5, 50)

ggsave(filename = "Total_activity_counts.pdf", plot = totalctplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Total Activity Duration ####
totaldurplot <- FinalData %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgtotaldur = mean(totaldur)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgtotaldur) %>%
  ggplot(aes(x = time, y = avgtotaldur)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Total activity duration") + 
  xlim(5, 50)

ggsave(filename = "Total_activity_duration.pdf", plot = totaldurplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)

#### Total Activity Distance ####
totaldistplot <- FinalData %>%
  group_by(Chemical, `Dose(mg/L)`, time) %>%
  mutate(avgtotaldist = mean(totaldist)) %>%
  ungroup() %>%
  group_by(Chemical) %>%
  select(Chemical, Group, `Dose(mg/L)`, time, avgtotaldist) %>%
  ggplot(aes(x = time, y = avgtotaldist)) +
  geom_line(aes(color = as.factor(`Dose(mg/L)`))) +
  geom_rect(data = dark, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = 'black', alpha = 0.2) +
  facet_wrap(~Chemical, scales = "free") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  xlab("Time [min]") +
  ylab("Total activity duration") + 
  xlim(5, 50)

ggsave(filename = "Total_activity_distance.pdf", plot = totaldistplot, device = "pdf", 
       path = paste0(getwd(), "/Images"), width = 1920, height = 1080,
       units = "px", scale = 3
)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                 Benchmark Dose                #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
