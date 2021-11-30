####Libraries####
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test
library(purrr)
library(ggfortify)

####Import metadata####
metadata <- read.csv("metadata.csv")
colnames(metadata) <- c("Chemical","10-1","10-2","10-3","10-4","10-5","3uM","Control")
chemnames <- metadata[1:6,1]
doses <- colnames(metadata[,2:8])

####Grouping####
##Chemical Groups
chemgroups <- c()
for(i in 1:length(chemnames)) {
  temp1 <- rowSums(metadata[i,2:8])
  temp2 <- rep(chemnames[i], temp1)
  chemgroups <- c(chemgroups, temp2)
  remove(temp1, temp2)
}

##Dose Groups
dosegroups <- c()
for(i in 1:length(chemnames)) {
  temp1 <- as.numeric(metadata[i,2:8])
  temp2 <- rep(doses, temp1)
  dosegroups <- c(dosegroups, temp2)
  remove(temp1, temp2)
}

##Dose and Chemical Groups Combined
dosechemgroups <- paste(chemgroups, dosegroups, sep = "-")

####Data Import####
raw_data <- list ()
for(i in chemnames) {
  filenames <- list.files(paste0(i))
  for(j in filenames){
    raw_data[[i]][[j]] <- read.table(paste0(i,"/",j),
                                header = FALSE,
                                stringsAsFactors = FALSE,
                                sep = "\t",
                                strip.white = TRUE) [-c(1:4),-c(2:3)]
    colnames(raw_data[[i]][[j]]) <- c("gene", j)
  }
}

####Adjust row and col names####
##Set gene col as row names
for(i in 1:length(raw_data)) {
  raw_data[[i]] <- raw_data[[i]] %>%
    reduce(left_join, by = "gene") %>%
    remove_rownames() %>%
    column_to_rownames("gene")
}
##Clean up col names
for(i in 1:length(raw_data)) {
  temp <- gsub(pattern = "-GQReadsPerGene.out.tab", 
               replacement = "", 
               x = names(raw_data[[i]]), 
               fixed = TRUE)
  colnames(raw_data[[i]]) <- temp
  remove(temp)
}

####Median Filtering####
##set threshold - For no median filtering, median_thrs <- 0
median_thrs <- 5
##Remove genes that do not meet median threshold
raw_data_median_filtered <- list()
for(i in chemnames) {
  raw_data_median_filtered[[i]] <- raw_data[[i]] %>%
    mutate(count_median = apply(raw_data[[i]], 1, median)) %>%
    filter(count_median >= median_thrs) %>%
    select(-count_median)
}

####Combined all the data together####
for(i in 1:length(raw_data_median_filtered)) {
  raw_data_median_filtered[[i]] <- raw_data_median_filtered[[i]] %>%
    rownames_to_column("gene")
}

my_merge <- function(df1, df2) {
  merge(df1, df2, by = "gene", all = TRUE)
}

combined_data_median_filtered <- Reduce(my_merge, raw_data_median_filtered) %>%
  remove_rownames() %>%
  column_to_rownames("gene")
combined_data_median_filtered[is.na(combined_data_median_filtered)] <- 0

####edgeR TMM Norm####
testDGE <- DGEList(counts = combined_data_median_filtered,
                   group = dosechemgroups,
                   genes = as.data.frame(row.names(combined_data_median_filtered)) %>% 
                     rename(genes = "row.names(combined_data_median_filtered)")
                   )
####GROUP NOT MATCHING...?
testDGE2 <- calcNormFactors(testDGE)
testDGE2$samples
testDGEtoDF <- edgeR::as.data.frame.DGEList(testDGE2)
##Not normalizing...?


####TESTING####
# test <- gsub(pattern = '_.*', replacement = "", x = colnames(combined_data_median_filtered))
# colnames(combined_data_median_filtered) <- test


####Filtering####
# mean_thrs <- 0
# median_thrs <- 0
# zero_prop_thrs <- 0.7 #proportion of samples with counts > 0
# #var_thrs <-  0.3
# #abundance_thrs <- 0.15 #filter out lowest proportion of genes (0.2 = 20%)



# filtered_data <- count_data %>%
#   mutate(count_mean = rowMeans(count_data)) %>%
#   mutate(count_median = apply(count_data, 1, median)) %>%
#   mutate(zero_proportion = rowSums(count_data != 0) / ncol(count_data)) %>%
#   mutate(count_total = rowSums(count_data)) %>%
#   mutate(count_vars = apply(count_data, 1, var)) %>%
#   filter(
#     count_mean >= mean_thrs,
#     count_median >= median_thrs,
#     zero_proportion >= zero_prop_thrs
#     ) %>%
#   #slice_max(order_by = count_vars, prop = 1 - var_thrs) %>%
#   #slice_max(order_by = count_total, prop = 1 - abundance_thrs) %>%
#   select(-count_median,
#          -count_mean,
#          -count_total,
#          -count_vars,
#          -zero_proportion) %>%
#   print()


####normalization####

####PCA####

test1 <- as.data.frame(t(log2(combined_data_median_filtered+1)))
chemgroupsdf <- as.data.frame(chemgroups)
test2 <- cbind(chemgroupsdf,test1)
#can't use chemgroups column... prcomp breaks

pca_model <- prcomp(test2)
autoplot(pca_model)

# library(FactoMineR)
# test.pca <- PCA(t(combined_data_median_filtered))
# plot(test.pca)

