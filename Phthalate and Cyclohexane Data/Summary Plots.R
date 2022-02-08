####Libraries####
library(tidyverse)
source("Bootstrapping_Functions.R")
source("mode_antimode.R")

####Data Import####

#Metadata
metadata <- read.csv("RNAseqData/metadata_nocontrol.csv")
chemnames <- unique(metadata$chemical)
lowestdoses <- unique(metadata[, c("chemical", "dose")]) %>%
  group_by(chemical) %>%
  filter(dose > 0) %>%
  summarise_all(min)

#RDS files
logBMDvalues <- readRDS("all_BMD_list_logtransformed.RDS")
goBMDvalues <- readRDS("go_BMD_list_logtransformed.RDS")
reactomeBMDvalues <- readRDS("reactome_BMD_list_logtransformed.RDS")

#Apical Data

apicallist <- read.csv("apical_endpoints.csv")

####Bootstrapping####
#20th gene
nthgenelist <- sapply(logBMDvalues, nth_gene_bootstrap) %>% 
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "nthgene") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

#10th percentile
nthpercentlist <- sapply(logBMDvalues, nth_percent_bootstrap) %>%
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "nthpercent") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

#first mode
modelist <- sapply(logBMDvalues, mode_bootstrap) %>%
  t() %>% 
  as.data.frame %>%
  mutate(endpoint = "mode") %>%
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
  mutate(endpoint = "goterm") %>%
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
  mutate(endpoint = "reactome") %>%
  rownames_to_column("chemical") %>%
  tibble() %>%
  relocate(chemical, endpoint)

remove(reactometemplist)




####Summary####



bigtibble <- bind_rows(nthgenelist, nthpercentlist, modelist, gotermlist, reactomelist) %>%
  rename(lowerCI = "2.5%", median = "50%", upperCI = "97.5%") %>%
  bind_rows(apicallist)

summaryplots <- list()

i <- "DBP"

# for(i in chemnames){
#   summaryplots[[i]] <-
    bigtibble %>%
    filter(chemical == i) %>%
    ggplot(aes(x=median, y=endpoint)) +
    xlim(-5,1) + 
    geom_point(size = 3, aes(colour = factor(endpoint)), show.legend = FALSE) +
    geom_errorbarh(aes(xmin=lowerCI, xmax=upperCI, colour = factor(endpoint)), height = .2, size = 0.5, show.legend = FALSE) +
    # geom_pointrange(aes(xmin = lowerCI, xmax = upperCI), size=0.5) +
    labs(title = paste0(i," (n=",length(logBMDvalues[[i]][[1]]), " genes)"), x = "logBMD", y = "Endpoint") +
    theme(plot.title = element_text(hjust = 0.5, size=12)) +
    theme_classic()
# }

source("BMDExpressFunctions.R")

multiplot(plotlist = summaryplots, cols = 3 )