####Libraries####
library(here)
library(tidyverse)
library(rlang)
library(Rcurvep)
source("cal_auc_simi_endpoints.R")

####Metadata and Import####
metadata <- read.csv("meta_data_behaviour.csv")
filenames <- metadata[,1]
chemnames <- substr(filenames, 1, nchar(filenames)-4)
HighDose <- setNames(metadata[,2], chemnames)
LowDose <- setNames(metadata[,3],chemnames)
raw_data <- list()
for(i in chemnames){
  raw_data[[i]] <- read.table(paste0(i, "/", i, ".txt"),
                              header = TRUE)
}

####Calculate Similarity Endpoints####
simi_endps <- list("ld_pearson" = seq(1, 20, by = 1))
simi_endpoints <- list()
for(i in chemnames){
  temp <- create_simi_endpoints(raw_data[[i]], segments = simi_endps, metric = "pearson")
  simi_endpoints[[i]] <- temp[[1]]
  rm(temp)
}

####Normalize Data####
#simi_normalize function - to be moved later
simi_normalize <- function(x){
  x %>%
    mutate(control_median = unname(tapply(x$endpoint_value, x$is_VC, median)["TRUE"])) %>%
    mutate(endpoint_value_norm = (endpoint_value/control_median)*100-100) %>%
    select(-c(endpoint_value, control_median)) %>%
    mutate(dose = substr(embryo_id, 1, nchar(embryo_id)-2)) %>%
    select(-embryo_id)
    # mutate(embryo_id = substr(embryo_id, nchar(embryo_id), nchar(embryo_id)))
}

simi_norm <- lapply(simi_endpoints, simi_normalize)

#dose_replacement function - to be moved later
dose_replacement <- function(x, Highdose, foldchange = 10){
  ndosegroups <- 1:(length(unique(x$dose))-1)
  for(i in ndosegroups){
    x$dose[x$dose == paste0("Dose", i)] <- Highdose/foldchange^(i-1)
  }
  x$dose[x$dose == "Control"] <- 0
  x$dose <- as.numeric(x$dose)
  return(x)
}

for(i in chemnames){
  simi_norm[[i]] <- dose_replacement(x = simi_norm[[i]], Highdose = HighDose[i]) %>%
    group_by(dose)
}

####ANOVA####
summarystats <- function(x, grouping = "dose", values = "endpoint_value_norm"){
  group_by(x, dose) %>%
    summarise(
      count = n(),
      mean = mean(endpoint_value_norm, na.rm = TRUE),
      sd = sd(endpoint_value_norm, na.rm = TRUE)
    )
}

summarystats_list <- lapply(simi_norm, summarystats)

combinedanova <- function(x){
  aov(endpoint_value_norm ~ dose, data = x) %>% 
    summary()
}

anova_list <- lapply(simi_norm, combinedanova)

for(i in chemnames){
  write.table(simi_norm[[i]], file = paste0(i,"/",i,"_simi_norm.txt"))
}



####Plotting####
behaviourplots <- list()

for(i in chemnames) {
  behaviourplots[[i]] <-
    simi_norm[[i]] %>%
    ggplot(aes(x = as.character(dose), y = endpoint_value_norm, group = dose)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.2,
                height = 0,
                colour = "black") +
    labs(title = paste0(i), x = "Dose (µg/L)", y = "Response") +
    theme_classic()
}

print(behaviourplots[[1]])

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(plotlist = behaviourplots, cols = 2 )









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

