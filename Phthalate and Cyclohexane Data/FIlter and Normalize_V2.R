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
chemdosegroups <- paste(chemgroups, dosegroups, sep = "_")

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

# ####Median Filtering####
# ##set threshold - For no median filtering, median_thrs <- 0
# median_thrs <- 5
# ##Remove genes that do not meet median threshold
# raw_data_median_filtered <- list()
# for(i in chemnames) {
#   raw_data_median_filtered[[i]] <- raw_data[[i]] %>%
#     mutate(count_median = apply(raw_data[[i]], 1, median)) %>%
#     filter(count_median >= median_thrs) %>%
#     select(-count_median)
# }

####Combined all the data together####
for(i in 1:length(raw_data)) {
  raw_data[[i]] <- raw_data[[i]] %>%
    rownames_to_column("gene")
}

my_merge <- function(df1, df2) {
  merge(df1, df2, by = "gene", all = TRUE)
}

combined_raw_data <- Reduce(my_merge, raw_data) %>%
  remove_rownames() %>%
  column_to_rownames("gene")

transposed_raw_data <- t(combined_raw_data)
grouping <- as.factor(chemdosegroups)
chemical <- as.factor(chemgroups)
dose <- as.factor(dosegroups)
temp <- data.frame(chemical, dose, transposed_raw_data) 
grouped_data <- temp %>% group_by(chemical, dose)
summary_data <- grouped_data %>% summarise_all(median)
summary_data_trimmed <- summary_data %>%
  data.frame() %>%
  select(-chemical,
         -dose) %>%
  summarise_all(min)
rm(temp)

medianthrs <- 5

unfiltered_data <- t(rbind(min_median = summary_data_trimmed, transposed_raw_data))
median_filtered_data <- unfiltered_data %>%
  data.frame() %>%
  filter(min_median >= medianthrs) %>%
  select(-min_median)


####Combined Normalization####
#Not very well done...
only_DMSO <- median_filtered_data[,16:19]
no_DMSO <- median_filtered_data[,c(1:16, 21:34, 39:52, 57:71, 76:92, 97:110)]
median_filtered_data_combined_control <- cbind(only_DMSO, no_DMSO)

newdoses <- colnames(metadata[,2:7])
#new dose groups
newdosegroups <- c()
for(i in 1:length(chemnames)) {
  temp1 <- as.numeric(metadata[i,2:7])
  temp2 <- rep(newdoses, temp1)
  newdosegroups <- c(newdosegroups, temp2)
  remove(temp1, temp2)
}
newdosegroups <- c(rep("control", 4), newdosegroups)

##Chemical Groups
newchemgroups <- c()
for(i in 1:length(chemnames)) {
  temp1 <- rowSums(metadata[i,2:7])
  temp2 <- rep(chemnames[i], temp1)
  newchemgroups <- c(newchemgroups, temp2)
  remove(temp1, temp2)
}
newchemgroups <- c(rep("DMSO", 4), newchemgroups)
newchemdosegroups <- paste(newchemgroups, newdosegroups, sep = "_")

sizeFactor <- calcNormFactors(median_filtered_data_combined_control)
norm_combined_data <- median_filtered_data_combined_control
for(i in 1:ncol(norm_combined_data)){
  norm_combined_data[,i]<-norm_combined_data[,i]/(sum(norm_combined_data[,i])*sizeFactor[i])
}
colnames(norm_combined_data) <- newchemdosegroups

####Separate Normalization####

chemfactor <- as.factor(chemgroups)
temp1 <- data.frame(chemfactor, t(median_filtered_data))
for(i in chemnames){
  assign(i, subset(temp1, chemfactor == i) %>% select(-chemfactor) %>% t())
}
rm(temp1)
for(i in chemnames){
  assign(paste("sizeFactor", i, sep = "_"), calcNormFactors(get(i)))
}
for(i in chemnames){
  assign(paste("norm_data", i, sep = "_"), data.frame(get(i)))
}

# for(i in chemnames){
#   for(j in paste("sizeFactor",i,sep = "_")){
#     for(k in paste("norm_data",i,sep = "_")){
#       for(l in 1:ncol(get(k))){
#         k[,l] <- k[,l]/(sum(k[,l])*j[l])
#       }
#     }
#   }
# } #gives back "Error in k[, l] : incorrect number of dimensions"

for(i in 1:ncol(norm_data_DBC)){
  norm_data_DBC[,i]<-norm_data_DBC[,i]/(sum(norm_data_DBC[,i])*sizeFactor_DBC[i])
}
for(i in 1:ncol(norm_data_DBP)){
  norm_data_DBP[,i]<-norm_data_DBP[,i]/(sum(norm_data_DBP[,i])*sizeFactor_DBP[i])
}
for(i in 1:ncol(norm_data_DEHC)){
  norm_data_DEHC[,i]<-norm_data_DEHC[,i]/(sum(norm_data_DEHC[,i])*sizeFactor_DEHC[i])
}
for(i in 1:ncol(norm_data_DEHP)){
  norm_data_DEHP[,i]<-norm_data_DEHP[,i]/(sum(norm_data_DEHP[,i])*sizeFactor_DEHP[i])
}
for(i in 1:ncol(norm_data_P1400)){
  norm_data_P1400[,i]<-norm_data_P1400[,i]/(sum(norm_data_P1400[,i])*sizeFactor_P1400[i])
}
for(i in 1:ncol(norm_data_SANT)){
  norm_data_SANT[,i]<-norm_data_SANT[,i]/(sum(norm_data_SANT[,i])*sizeFactor_SANT[i])
}



##testing
# raw_data2 <- lapply(raw_data, column_to_rownames, "gene")
# list2env(raw_data2, envir = .GlobalEnv)
# for(i in raw_data2){
#   temp1 <- calcNormFactors(i)
#   temp2 <- i
#   for(j in 1:ncol(temp2)){
#     temp2[,j] <- temp2[,j]/sum(temp2[,j]*temp1[j])
#   }
# }

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

temp1 <- as.data.frame(t(norm_combined_data))
chemgroupsfactor <- as.factor(newchemgroups)
grouped_norm_data <- data.frame(chemgroupsfactor,temp1)
rm(temp1)
pca_model <- prcomp(grouped_norm_data[,-1])
autoplot(pca_model, colour = 'chemgroupsfactor', data = grouped_norm_data, label = TRUE)




