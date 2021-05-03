setwd("E:\\Behaviour\\Behaviour R Outputs")
library(tidyr)
library(tidyverse)

#load file
data1 <- read.table("E:\\Behaviour\\TGSH\\TGSH.XLS", header = TRUE)

data1 <- data1 %>%
  drop_na() %>% #drop all rows with NA
  filter(row_number() %% 2 != 0) %>% #filter and remove every even row. For x != y, x= how often to delete and y= where to begin.
  subset(, select = -c(animal, an)) #remove "animal" and "an" column

