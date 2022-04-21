####Libraries####
library(tidyverse)
library(biomaRt)
source("BMDExpressFunctions.R")
options(scipen = 9)

####Metadata Import####
metadata <- read.csv("RNAseqData/metadata.csv")
chemnames <- unique(metadata$chemical)

lowestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(min)

highestdoses <- unique(metadata[,c("chemical","dose")]) %>%
  group_by(chemical) %>%
  filter(dose>0) %>%
  summarise_all(max)

####Data Import####
raw_BMD<- list()
filenames <- list.files("BMDExpressData/BMD")

for(i in 1:length(chemnames)){
  raw_BMD[[chemnames[i]]] <- read.table(paste0("BMDExpressData/BMD/",filenames[i]), header = TRUE, sep = "\t") %>%
    cleanupcolumns()
}

# Filter
filtered_BMD <- list()
for(i in chemnames){
  filtered_BMD[[i]] <- BMDfiltering(x = raw_BMD[[i]], 
                                          lowdose = lowestdoses$dose[lowestdoses$chemical==i], 
                                          highdose = highestdoses$dose[highestdoses$chemical==i],
  ) %>%
    mutate(logBMD = log10(BMD))
}


####Gene List Import####
geneset <- read.csv("Custom Gene Sets/estrogen_genes.csv", stringsAsFactors = F) %>%
  pull(gene)

#Biomart#
#import Biomart dataset
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('drerio_gene_ensembl', mart)

#retrieve gene symbol and ensembl id
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'zfin_id_symbol',
    'ensembl_gene_id'),
  uniqueRows = TRUE)

#Filter gene set#
ensembl_geneset <- annotLookup %>%
  filter(zfin_id_symbol %in% geneset)

####Filter BMD based on Gene Set####
geneset_BMD <- list()

for(i in chemnames){
  temp <- filtered_BMD[[i]] %>%
    rename(ensembl_gene_id = Probe.Id)
  geneset_BMD[[i]] <- inner_join(ensembl_geneset, temp,  by = "ensembl_gene_id") %>%
    rename("ZFIN Gene ID" = zfin_id_symbol,
           "Ensembl Gene ID" = ensembl_gene_id) %>%
    dplyr::select(-logBMD,
                  -fitPValue)
  rm(temp)
}

geneset_BMD_df <- bind_rows(geneset_BMD, .id = 'Chemical')

write.csv(geneset_BMD_df, "BMDExpressData/Output/geneset_BMD.csv", row.names = F)
