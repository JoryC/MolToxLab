#Library
library(tidyverse)
library(DT)
source("BMDExpressFunctions.R")

#Import DEGs
files <-
  list.files(path = "~/MolToxBackups/Jory WT 10 Chemicals/RNAseqData/output", pattern = "DEG_table.txt")
DEGs_list <- NULL
for (i in 1:length(files)) {
  temp_file <- files[i]
  temp <-
    read.table(
      file = paste0(getwd(), "/RNAseqData/output/", files[i]),
      sep = "\t",
      header = TRUE
    )
  temp <- temp %>%
    rownames_to_column(var = "gene")
  DEGs_list[[str_split(temp_file, pattern = "_", simplify = TRUE)[, 1]]] <-
    temp
}

#Plot DEG counts per chem
p_list <- NULL
n_DEGs_list <- NULL
for (i in 1:length(DEGs_list)) {
  n_DEGs <- DEGs_list[[i]] %>%
    mutate(gene = str_split(gene, pattern = "[.]", simplify = TRUE)[, 1]) %>%
    select(gene) %>%
    unique() %>%
    count()
  p <- DEGs_list[[i]] %>%
    mutate(comparison = str_split(
      string = comparison,
      pattern = "_",
      simplify = TRUE
    )[, 2]) %>%
    mutate(comparison = factor(
      comparison,
      levels = unique(comparison),
      ordered = TRUE
    )) %>%
    group_by(comparison) %>%
    unique() %>%
    tally(name = "DEG_count") %>%
    ggplot(aes(x = comparison, y = DEG_count)) +
    geom_bar(stat = "identity", fill = "black") +
    #geom_smooth(aes(group = 1), method = "lm", se= FALSE, color = 'red', size = 1) +
    theme_classic() +
    xlab("Dose (mg/L)") +
    ylab("Number of DEGs") +
    ggtitle(names(DEGs_list[i]), subtitle = paste0("n=", n_DEGs, " genes"))
  
  
  p_list[[names(DEGs_list[i])]] <- p
  n_DEGs_list[[names(DEGs_list[i])]] <- n_DEGs
}

#Plot
multiplot(plotlist = p_list,
          layout = matrix(
            c(1:length(p_list), NA, NA),
            nrow = 3,
            ncol = 3,
            byrow = TRUE
          ))

#Pull the DEG tables for supplementary Information
all_DEGs <- plyr::ldply(DEGs_list) %>%
  mutate(gene = str_split(gene, pattern = "[.]", simplify = TRUE)[, 1]) %>%
  mutate(comparison = str_split(comparison, pattern = "_", simplify = TRUE)[, 2])
datatable(
  all_DEGs,
  rownames = FALSE,
  fillContainer = TRUE,
  extensions = 'Buttons',
  options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    lengthMenu = list(c(10, 25, 50, -1),
                      c(10, 25, 50, "All"))
  )
)
