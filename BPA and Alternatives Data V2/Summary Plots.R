####Libraries####
library(tidyverse)
source("Bootstrapping_Functions.R")
source("mode_antimode.R")
source("BMDExpressFunctions.R")
options(scipen = 9)

####Data Import####

#Metadata
metadata <- read.csv("RNAseqData/metadata.csv")
chemnames <- unique(metadata$chemical)

lowestdoses <- unique(metadata[, c("chemical", "dose")]) %>%
  group_by(chemical) %>%
  filter(dose > 0) %>%
  summarise_all(min)

highestdoses <- unique(metadata[, c("chemical", "dose")]) %>%
  group_by(chemical) %>%
  filter(dose > 0) %>%
  summarise_all(max)

#RDS files
logBMDvalues <- readRDS("BMDExpressData/RDS/all_BMD_list_logtransformed.RDS")
goBMDvalues <- readRDS("BMDExpressData/RDS/go_BMD_list_logtransformed.RDS")
reactomeBMDvalues <- readRDS("BMDExpressData/RDS/reactome_BMD_list_logtransformed.RDS")

#Apical Data

apicallist <- read.csv("apical_endpoints.csv")

####Bootstrapping####
#20th gene
nthgenelist <- sapply(logBMDvalues, nth_gene_bootstrap) %>% 
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "20th gene") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

#10th percentile
nthpercentlist <- sapply(logBMDvalues, nth_percent_bootstrap) %>%
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "10th pct") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

#first mode
modelist <- sapply(logBMDvalues, mode_bootstrap) %>%
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "1st Mode") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)


####GO Term####

gotemplist <- list()
for(i in chemnames){
  for(j in 1:length(goBMDvalues[[i]])){
    if(length(goBMDvalues[[i]]) > 0){
      gotemplist[[i]][[j]] <- pathway_bootstrap(goBMDvalues[[i]][j])
    } else {
      gotemplist[[i]] <- list
      print(paste("Missing values in", i))
    }
  }
}

gotermlist <- sapply(gotemplist, averageCI) %>%  
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "GO Term") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

remove(gotemplist)

####Reactome####

reactometemplist <- list()
for(i in chemnames){
  for(j in 1:length(reactomeBMDvalues[[i]])){
    if(length(reactomeBMDvalues[[i]]) > 0){
      reactometemplist[[i]][[j]] <- pathway_bootstrap(reactomeBMDvalues[[i]][j])
    } else {
      reactometemplist[[i]] <- list()
      print(paste("Missing values in", i))
    }
  }
}

reactomelist <- sapply(reactometemplist, averageCI) %>%
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "Reactome") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

remove(reactometemplist)

####Summary####

bigtibble <- bind_rows(nthgenelist, nthpercentlist, modelist, gotermlist, reactomelist) %>%
  rename(lowerCI = "2.5%", median = "50%", upperCI = "97.5%") %>%
  bind_rows(apicallist) %>% 
  mutate(missingvalues = if_else(is.na(median), "ND", "")) %>%
  mutate(endpoint = factor(endpoint, levels = c("Reactome",
                                                "GO Term",
                                                "10th pct",
                                                "20th gene",
                                                "1st Mode",
                                                "Behaviour",
                                                "aPOD")))

summaryplots <- list()

# i <- "DBP"

for(i in chemnames){
  summaryplots[[i]] <-
    bigtibble %>%
    filter(chemical == i) %>%
    ggplot(aes(x=median, y=endpoint)) +
    xlim(-5,3) +
    # scale_x_continuous(limit = c(-5, 2), oob = function(x, limits) x) +
    geom_point(size = 3, aes(colour = endpoint), show.legend = FALSE) +
    geom_errorbarh(aes(xmin=lowerCI, xmax=upperCI, colour = endpoint), height = .3, size = 0.5, show.legend = FALSE) +
    # geom_pointrange(aes(xmin = lowerCI, xmax = upperCI), size=0.5) +
    geom_text(aes(x = -1, label = missingvalues)) +
    labs(title = paste0(i," (n=",length(logBMDvalues[[i]][[1]]), " genes)"), x = "logBMD (µM)", y = "Endpoint") +
    theme(plot.title = element_text(hjust = 0.5, size=12)) +
    theme_classic()
}

multiplot(plotlist = summaryplots, cols = 3 )