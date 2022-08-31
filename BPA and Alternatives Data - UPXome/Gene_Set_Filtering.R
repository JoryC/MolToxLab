####Libraries####
library(biomaRt)
library(tidyverse)
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

####Fisher's Exact Test####

#Import normalized Data
normdata_gene <- list()

for(i in chemnames){
  normdata_gene[[i]] <- read.table(paste0("RNAseqData/normalizedData/", i, "_normData.txt"),
                                   header = TRUE) %>%
  select(gene)
}

#Total number of estrogen genes in normalized Data
geneset_in_norm <- list()

for(i in chemnames){
  geneset_in_norm[[i]] <- normdata_gene[[i]] %>%
    filter(gene %in% (ensembl_geneset %>% pull(ensembl_gene_id)))
}

#Total number of estrogen genes without BMD
geneset_noBMD <- list()

for(i in chemnames){
  geneset_noBMD[[i]] <- geneset_in_norm[[i]] %>%
    filter(!gene %in% (geneset_BMD[[i]] %>% pull("Ensembl Gene ID")))
}

#Total number of non-estrogen genes with BMD
non_geneset_BMD <- list()

for(i in chemnames){
  non_geneset_BMD[[i]] <- filtered_BMD[[i]] %>%
    select(Probe.Id) %>%
    filter(!Probe.Id %in% (geneset_BMD[[i]] %>% pull("Ensembl Gene ID")))
}

#Total number of non-estrogen genes without BMD
non_geneset_noBMD <- list()

for(i in chemnames){
  non_geneset_noBMD[[i]] <- normdata_gene[[i]] %>%
    filter(!gene %in% (ensembl_geneset %>% pull(ensembl_gene_id))) %>%
    filter(!gene %in% (filtered_BMD[[i]] %>% pull("Probe.Id")))
}

#Fisher's Test Table Setup
fisher_table <- list()

for(i in chemnames){
  fisher_table[[i]] <- data.frame(
    "GeneSet" = c(nrow(geneset_BMD[[i]]), nrow(geneset_noBMD[[i]])),
    "Non_Geneset" = c(nrow(non_geneset_BMD[[i]]), nrow(non_geneset_noBMD[[i]])),
    row.names = c("BMD", "No BMD"),
    stringsAsFactors = F
  )
}

#Fisher's Test and p value
fisher_test_list <- list()
fisher_p_value <- list()

for(i in chemnames){
  fisher_test_list[[i]] <- fisher.test(fisher_table[[i]])
  fisher_p_value[[i]] <- fisher_test_list[[i]]["p.value"] %>% 
    as.numeric()
}


#Geneset BMD
geneset_BMD_df <- bind_rows(geneset_BMD, .id = 'Chemical')

for(i in 1:nrow(geneset_BMD_df)){
  geneset_BMD_df[i,"p-value"] <- fisher_p_value[geneset_BMD_df[i,"Chemical"]]
}

geneset_BMD_df <- geneset_BMD_df %>%
  group_by(Chemical) %>%
  summarise(across("ZFIN Gene ID", ~paste(., collapse = ", ")),
            across(c("BMD", "BMDL", "BMDU"), median))

write.csv(geneset_BMD_df, "BMDExpressData/Output/geneset_BMD.csv", row.names = F)
