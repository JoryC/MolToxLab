#Load the required packages
library(readr)
library(dplyr)
library(purrr)

#### BPAF ####

#Raw data has doses in row format. We want to switch that to columns using the t() function in our call to read_delim()
BPAF_Base <- as.data.frame(t(read_delim(file = "~/MolToxLab/Alamar_Blue/BPAF_2021_03_01/BPAF_Baseline_Cut.csv", #Input your file path and state your object name
                                    delim = "\t", #tab separated
                                    col_names = as.character(1:9) #imports 9 rows of your cut file. Please see shell script txt_2_csv for cutting txt files
      )
    )
  )
#Assign Column names to sort later. The second argument in rep() tells us how many times we ant to replicate the string of characters. We have 2 plates for this chemical so we will use 2
colnames(BPAF_Base) <- rep(c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Carrier Control (0.01%DMSO)", "Empty Control"), 2)

View(BPAF_Base)

#Combine all of the duplicate columns into a new tibble
BPAF_Base <- BPAF_Base %>%
  map(as.character) %>% #creates a list of column names
  split(names(.)) %>% #organizes names
  map_df(flatten_chr) %>% #unlist and assemble into df
  mutate_all(function(x) as.numeric(as.character(x)))
glimpse(BPAF_Base) #verify that your values are <dbl>

#Clean up the Empty Control Wells
BPAF_Base[c(5:9,14:18),7] = NA
View(BPAF_Base)

#Save a copy of your work
write.table(BPAF_Base, file = "~/MolToxLab/Alamar_Blue/BPAF_2021_03_01//BPAF_Base.txt", sep = "\t", row.names = FALSE)

#### BPA ####

#Raw data has doses in row format. We want to switch that to columns using the t() function in our call to read_delim()
BPA24 <- as.data.frame(t(read_delim(file = "~/MolToxLab/Alamar_Blue/BPA_2021_02_15/BPA24h_Cut.csv", #Input your file path
                                     delim = "\t", #tab separated
                                     col_names = as.character(1:9) #imports 9 rows
    )
  )
)
#Assign Column names to sort later. The second argument in rep() tells us how many times we ant to replicate the string of characters. We have 3 plates for this chemical so we will use 3
colnames(BPA24) <- rep(c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Carrier Control (0.01%DMSO)", "Empty Control"), 3)

View(BPA24)

#Combine all of the duplicate columns into a new tibble
BPA24 <- BPA24 %>%
  map(as.character) %>% #creates a list of column names
  split(names(.)) %>% #organizes names
  map_df(flatten_chr) %>% #unlist and assemble into df
  mutate_all(function(x) as.numeric(as.character(x)))
glimpse(BPA24) #verify that your values are <dbl>

#Save a copy of your work
write.table(BPA24, file = "~/MolToxLab/Alamar_Blue/BPA_2021_02_15/BPA24.txt", sep = "\t", row.names = FALSE)

#Clean up the Empty Control Wells
BPA24[c(5:9,14:18,23:27),7] = NA
View(BPA24)
