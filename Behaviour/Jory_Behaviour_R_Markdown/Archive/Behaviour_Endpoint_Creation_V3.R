####Libraries####
library(here)
library(tidyverse)
library(rlang)
library(Rcurvep)
library(DescTools)
library(data.table)
library(plyr)
source("Functions/cal_auc_simi_endpoints.R")
source("Functions/behavioural_endpoint_calc.R")
options(scipen = 9)

####Metadata and Import####
metadata <- read.csv("Data/meta_data_behaviour_all.csv") #Metadata contains filenames, highest dose and lowst dose
chemnames <- list.files("Data/Sub") #list all the folders in the sub-directory inside /Data
filenames <- list.files(path = paste0("Data/Sub/", chemnames, "/"), pattern = "*.csv") #List the filenames in each subfolder
HighDose <- setNames(metadata[,2], chemnames) #subsetting metadata


#Import the raw unprocesed data
raw_data <- list()
for (i in 1:length(filenames)) {
  raw_data[[i]] <- read_delim(file = paste0("Data/All/", filenames[[i]]),
                          delim = "\t",
                          col_names = TRUE,
                          col_types = "cclnninninninnin", #Where c is character, l is logical, n is numeric, i is integer
                          col_select = -c("inadist", "emptyct", "emptydur"), #get rid of meaningless variables
                          na = c("", 'NA', "NA", "\t NA", "\tNA") #Specify what to consider NA
  )
}
names(raw_data) <- chemnames


#Start processing
# Tidying and cleaning the messy raw data output from the Zebrabox
temp <- lapply(raw_data, FUN = drop_na) #Getting rid of empty wells

#Create one big data frame and add the Chemical name as a variable
temp <- rbindlist(temp, idcol = "plate_id")
temp <- temp %>%
  filter(an == TRUE, end <= 3000) %>%
  mutate(Treatment = rep(
    c(
      "Dose1_A","Dose1_B","Dose1_C","Dose2_A","Dose2_B","Dose2_C",
      "Dose3_A","Dose3_B","Dose3_C","Dose4_A","Dose4_B","Dose4_C",
      "Dose5_A","Dose5_B","Dose5_C","Control_A","Control_B","Control_C"
    ),
    times = 1300,
    each = 3
  )) %>% #Separating Dose from Group letter with an underscore
  separate(col = Treatment,
           into = c('Dose', 'Group'),
           convert = TRUE) %>% #splitting up the variables
    mutate(
      is_VC = if_else(
        condition = Dose == "Control",
        true = 1,
        false = 0
      ),
      totaldist = smldist + lardist,
      totaldur = smldur + lardur,
      totalct = smlct + larct,
      time_end = end / 60,
      velocity = totaldist / time_end,
      value = totaldist, #Choose from totaldist, totaldur, totalct, velocity, etc... variables in the raw_data dataframe... this is the endpoint we are analyzing
      animal = as.numeric(str_extract_all(animal, "[0-9]+")),
      embryo_id = paste(Dose, animal, sep = "_")
    ) %>%
  filter(time_end >= 21) %>%
  mutate(time_end = time_end-20) %>%
  select(plate_id, embryo_id, is_VC, time_end, value)


raw_data <- split(as.data.frame(temp), ~ plate_id) #Converting to a nested tibble for the rest of the pipeline


####Calculate Similarity Endpoints####
simi_endps <- list("ld_pearson" = seq(1, 30, by = 1))
simi_endpoints <- list()
for(i in chemnames){
  temp <- create_simi_endpoints(raw_data[[i]], segments = simi_endps, metric = "pearson")
  simi_endpoints[[i]] <- temp[[1]]
  rm(temp)
}


####Normalize Data####
simi_norm <- lapply(simi_endpoints, simi_normalize)

for(i in chemnames){
  simi_norm[[i]] <- dose_replacement(x = simi_norm[[i]], Highdose = HighDose[i]) %>%
    group_by(dose)
}

####ANOVA####
for(i in chemnames){
  simi_norm[[i]]$dose <- as.factor(simi_norm[[i]]$dose)
}

summarystats_list <- lapply(simi_norm, summarystats)

anova_list <- sapply(simi_norm, combinedanova)
print(anova_list)

dunnett_list <- sapply(simi_norm, combineddunnett) %>%
  setNames(., chemnames) %>% 
  as.array()
print(dunnett_list)

####Export Data####
for(i in chemnames){
  write.csv(simi_norm[[i]], file = paste0("Output/", i,"/",i,"_simi_norm.csv"))
}
for(i in chemnames){
  write.csv(anova_list[[i]], file = paste0("Output/", i,"/",i,"_anova.csv"))
}
for(i in chemnames){
  write.csv(dunnett_list[[i]], file = paste0("Output/", i,"/",i,"_dunnett.csv"))
}

simi_norm_tib <- as_tibble(rbindlist(simi_norm))
write_csv(simi_norm_tib, file = "Output/simi_norm_data.csv")

simi_norm_tib <- read_csv(file = "Output/simi_norm_data.csv")


###Find the Interesting Data###

#ANOVA results
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
#Yep, 34DCA is significant... be cautious look at the raw data

#write.csv(x = anova_comb, file = "Output/ANOVA_Results_Simi_endpoint.csv")

#Dunnett's Test rquires a bunch of wrangling
dunnett_temp <- list() #Take Dunnett's test results without any of the fancy summary information and shove it into a named list
for (i in chemnames) {
  dunnett_temp[[names(dunnett_list[i])]] <-
    dunnett_list[[i]] %>% as.data.frame() %>% rownames_to_column(var = "dose") #Coerce to a data frame temporarily so what we can take the row names of the reults and turn them into a variable with rownames_to_column
}
dunnett_comb <- ldply(dunnett_temp) #this function combines all of the lists together and gives them a variable name according to the chemical
dunnett_comb$dose = substr(dunnett_comb$dose,
                           start = 1,
                           stop = nchar(dunnett_comb$dose) - 2) %>%
  as.numeric() #Here we are fixing the dose column... the dose column has the test dose related to the control... but we just want to see what the test dose is without it giving us redundant information about the comparison to the control for every observation...
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
#Yep, dose 1000ug/L in 34DCA once again... but the raw data is no good in thehighest dose :(

#write.csv(x = dunnett_comb, file = "Output/Dunnett_Results_Simi_endpoint.csv")

####Plotting####
behaviourplots <- list()
for(i in chemnames) {
  behaviourplots[[i]] <-
    simi_norm[[i]] %>%
    ggplot(aes(x = as.character(dose), y = endpoint_value_norm, group = dose)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(position = position_jitter(width = 0.2,
                height = 0, seed = 42069),
                colour = "black") +
    labs(title = paste0(i)) + # x = "Dose (µg/L)", y = "Response") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
}
print(behaviourplots[[8]])

multiplot(plotlist = behaviourplots, cols = 5 )

#Playing with plots
 behaviourplots <- list()
 for(i in chemnames) {
   behaviourplots[[i]] <-
     simi_norm[[i]] %>%
     ggplot(aes(x = as.numeric(dose), y = endpoint_value_norm)) +
     geom_boxplot(aes(group = dose), outlier.shape = NA, width = 0.5) +
     geom_smooth(se = FALSE, method = "auto") +
     geom_jitter(width = 0.2,
                 height = 0,
                 colour = "black") +
     labs(title = paste0(i), x = "Dose (µg/L)", y = "Response") +
     scale_x_continuous(breaks = c(1:6), labels = levels(simi_norm[[i]]$dose)) +
     theme_classic()
 }
 print(behaviourplots[[1]])
 
 #Plotting from tibble
 simi_norm_tib %>%
   mutate(dose_num = as.numeric(as.vector(dose))) %>%
   group_by(plate_id) %>%
   ggplot(aes(x = dose_num, y = endpoint_value_norm, group = dose_num)) +
   geom_boxplot(outlier.shape = NA, width = 0.5) +
   geom_jitter(position = position_jitter(width = 0.2,
                                          height = 0, seed = 42069),
               colour = "black") +
   labs(x = "Dose (µg/L)", y = "Response") +
   theme_classic() +
   scale_x_log10() +
   facet_wrap(~plate_id)

####Setup Data for RCurvep####
#pre_rcurvep function - to be moved later
pre_rcurvep <- function(x){
  x %>%
    rename(
      resp = endpoint_value_norm,
      chemical = plate_id,
    ) %>%
    filter(is_VC == 0) %>%
    mutate(conc = dose) %>%
    select(-dose, -is_VC)
}

simi_prercurvep <- lapply(simi_norm, ungroup)
simi_prercurvep <- lapply(simi_prercurvep, pre_rcurvep)
simi_combined <- bind_rows(simi_prercurvep)

test <- estimate_dataset_bmr(combi_run_rcurvep(
  simi_combined, 
  n_samples = 10, 
  keep_sets = c("act_set"), 
  TRSH = seq(5, 95, by = 5) # test all candidates, 5 to 95
), plot = TRUE)

