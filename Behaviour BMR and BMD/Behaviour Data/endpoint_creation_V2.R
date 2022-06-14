####Libraries####
library(here)
library(tidyverse)
library(rlang)
library(Rcurvep)
library(DescTools)
source("Functions/cal_auc_simi_endpoints.R")
source("Functions/behavioural_endpoint_calc.R")
options(scipen = 9)

####Metadata and Import####
metadata <- read.csv("RawData/meta_data_behaviour.csv")
filenames <- metadata[,1]
chemnames <- substr(filenames, 1, nchar(filenames)-4)
HighDose <- setNames(metadata[,2], chemnames)
LowDose <- setNames(metadata[,3], chemnames)
raw_data <- list()
for(i in chemnames){
  raw_data[[i]] <- read.table(paste0("RawData/", i, "/", i, ".txt"),
                              header = TRUE)
}

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
    labs(title = paste0(i), x = "Dose (µg/L)", y = "Response") +
    theme_classic()
}

# print(behaviourplots[[1]])

multiplot(plotlist = behaviourplots, cols = 2 )

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

