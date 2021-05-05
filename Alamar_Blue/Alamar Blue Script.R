#### LOAD LIBRARIES ####
library(tidyverse)
library(abind)
library(car)

#### LOAD FILE AND CLEAN UP ####

#load files
baseline_1 <- read.table("Data\\TGSH_2021_02_02\\TGSH_Baseline_1.txt", skip = 13, nrows = 8, fill = TRUE)
baseline_2 <- read.table("Data\\TGSH_2021_02_02\\TGSH_Baseline_2.txt", skip = 13, nrows = 8, fill = TRUE)
after24h_1 <- read.table("Data\\TGSH_2021_02_02\\TGSH_24h_1.txt", skip = 13, nrows = 8, fill = TRUE)
after24h_2 <- read.table("Data\\TGSH_2021_02_02\\TGSH_24h_2.txt", skip = 13, nrows = 8, fill = TRUE)

#Remove last 3 columns and 7th row function
cleanup <- function(x) {
  x <- x[, -c(10:12)]
  x <- x[-7,]
  return(x)
}

#turn into a list to lapply cleanup function on all objects
templist <- list(baseline_1, baseline_2, after24h_1, after24h_2)
templist2 <- lapply(templist, cleanup)

#extract list into dataframes
baseline_1 <- data.frame(templist2[1])
baseline_2 <- data.frame(templist2[2])
after24h_1 <- data.frame(templist2[3])
after24h_2 <- data.frame(templist2[4])

rm(templist, templist2) #remove templists

#convert to array and get the mean of each cell
#along is along the "3rd dimension", i.e. not along rows or columns
baseline_avg <- rowMeans(abind(baseline_1, baseline_2, along = 3), dims=2)
after24h_avg <- rowMeans(abind(after24h_1, after24h_2, along = 3), dims=2)

#convert back to df
baseline_avg <- as.data.frame(baseline_avg)
after24h_avg <- as.data.frame(after24h_avg)

rm(after24h_1, after24h_2, baseline_1, baseline_2) #remove old separated data

#seperate out data from control (i.e. final row) NOTE:Breaks if run multiple times...
baseline_avg_control <- baseline_avg[7,1:4]
baseline_avg_control <- t(baseline_avg_control)
baseline_avg <- baseline_avg[-7,]
after24h_avg_control <- after24h_avg[7,1:4]
after24h_avg_control <- t(after24h_avg_control)
after24h_avg <- after24h_avg[-7,]

#assign highest dose and reformatting data #

#dose assignment
highestdose <- 1 #in mg/L
dose1 <- highestdose
dose2 <- highestdose/10
dose3 <- highestdose/100
dose4 <- highestdose/1000
dose5 <- highestdose/10000
dose6 <- 0

#formatting data into columns function
formatting <- function(x) {
  value <- as.vector(t(x)) #transpose data into a vector
  temp1 <- data.frame(value) #convert into a dataframe
  #add doses, rep function repeats string, "each" is how often each value is repeated before going to the next one
  temp1$dose <- rep(c(dose1, dose2, dose3, dose4, dose5, dose6), each = 9) 
  #add reps, rep function repeats string, "each" is how often each value is repeated before going to the next one
  temp1$replicate <- rep(c("A", "B", "C"), each = 1)
  return(temp1)
}

baseline <- formatting(baseline_avg)
after24h <- formatting(after24h_avg)
rm(after24h_avg, baseline_avg) #remove old dfs

#### NORMALIZE ####

#test data
deltavalue <- after24h$value - baseline$value
temp1 <- as.data.frame(cbind(deltavalue))
temp1$dose <- rep(c(dose1, dose2, dose3, dose4, dose5, dose6), each = 9) 
temp1$replicate <- rep(c("A", "B", "C"), each = 1)
rm(deltavalue)

#### ANCOVA ####

#set factor levels from highest to lowest
temp1$dose <- ordered(temp1$dose,
                         levels = c(dose1, dose2, dose3, dose4, dose5, dose6))

#summary of data
summary1 <- group_by(temp1, dose) %>% #(dataframe, group)
  summarise(
    count = n(),
    mean = mean(deltavalue, na.rm = TRUE), #replace deltavalue with data, na.rm is NA skipping
    sd = sd(deltavalue, na.rm = TRUE) #replace deltavalue with data, na.rm is NA skipping
  )
#check for homogeneity of variance
HomoVar <- leveneTest(deltavalue ~ dose, temp1) #should be insig, if not, not good

#run ancova, replicate is covariable, should not use type 1 since multiple factors
anovadose <- aov(deltavalue ~ dose + replicate, data = temp1) #aov(values ~ groups + covariable, data = x)
Anova(anovadose, type = "3") #intercept being significant just means grand mean =/= 0


#### PLOTTING ####
# (still testing)

temp2 <- temp1
temp2$mean <- rep((summary1$mean), each = 9)
temp2$sd <- rep((summary1$sd), each = 9)

ggplot() +
  geom_point(data = temp2, aes(dose, deltavalue)) +
  geom_point(data = temp2, aes(dose, mean), size = 3, color = 'red') +
  geom_errorbar(
    data = temp2,
    aes(dose, mean, ymin = mean - sd, ymax = mean + sd),
    width = 0.4,
    color = 'red'
  )
