library(dplyr)
library(ggplot2)
library(scales)
# library(plyr)
library(broom)
library(car)
library(DT)
library(outliers)
library(DescTools)
### FOLDERS, FILES AND FUNCTIONS
# setwd("~/MolToxLab/Lethality_and_Deformities/") # Please set your own working directory
temp1 <- list.files("Data/Apical_Data_Sheets/", "*.csv") #List all of the files in the working directory
apicaldata <- lapply(paste0("Data/Apical_Data_Sheets/", temp1), read.csv) #Import the data to a single object
rm(temp1)
chemnames <- gsub("\\.csv$", "", list.files("Data/Apical_Data_Sheets/", "*.csv"))
names(apicaldata) <- chemnames

#Fix Hatch Rate denominator to the # of alive fish instead of 20
lapply(X = apicaldata, FUN = function(x){mutate(x, Hatch.Rate = Hatched/Alive)})

#Output Sheets
all_apicaldata <- plyr::ldply(apicaldata) %>%
  rename(Chemical = .id)
write_csv(all_apicaldata, file = "master_apical_data.csv")


####Outliers####
grubbs.flag <- function(x) {
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  while(pv < 0.05) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(Affected.Rate=x,Outlier=(x %in% outliers)))
}

#Testing for Outliers...
input_apicaldata <- all_apicaldata %>%
  group_by(Chemical) %>%
  select(Affected.Rate) %>%
  mutate(Affected.Rate = round(Affected.Rate, digits = 3)) %>% #Bug with grubbs and large decimal places
  # filter(Chemical %in% c("BPA", "BPAF")) %>% #Testing with smaller data set
  nest()

grubbs_results <- input_apicaldata %>%
  mutate(grubbs_results = map(data, ~grubbs.flag(x = .x[[1]])))

apicaldata_join <- grubbs_results %>%
  select(-data) %>%
  unnest(cols = c(grubbs_results))

apicaldata_test <- all_apicaldata %>%
  mutate(Affected.Rate = round(Affected.Rate, digits = 3)) %>%
  inner_join(apicaldata_join, by = c("Chemical", "Affected.Rate"))
apicaldata_test <- apicaldata_test[-which(duplicated(apicaldata_test)),] #Delete duplicates introduced by join

test <- apicaldata_test %>%
  group_by(Chemical, Dose..mg.L.) %>%
  tally()

#Outliers removed
apicaldata_ol_rm <- apicaldata_test %>%
  mutate(across(
    .cols = c("Affected.Rate"),
    .fns = ~ replace(x = ., list = Outlier, values = NA) 
  ))


#nest/list data again
apicaldata <- apicaldata_ol_rm %>%
  group_by(Chemical) %>%
  split(f = .$Chemical)


####ANOVA####
#Survival Rate
anova_output_SR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Survival.Rate, data = apicaldata[[i]])
  anova_output_SR[i] <- summary(temp2)
  rm(temp2)
}
names(anova_output_SR) <- chemnames
print(anova_output_SR)
anova_SR_summary <- plyr::ldply(anova_output_SR, broom::tidy) %>%
  mutate(adj.p.value = p.adjust(p.value, "fdr"))
write_csv(anova_SR_summary, file = "anova_SR_summary.csv")

#Deformity Rate
anova_output_DR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Deformity.Rate, data = apicaldata[[i]])
  anova_output_DR[i] <- summary(temp2)
  rm(temp2)
}
names(anova_output_DR) <- chemnames
print(anova_output_DR)
anova_DR_summary <- plyr::ldply(anova_output_DR, broom::tidy) %>%
  mutate(adj.p.value = p.adjust(p.value, "fdr"))
write_csv(anova_DR_summary, file = "anova_DR_summary.csv")

#Hatch Rate
anova_output_HR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Hatch.Rate, data = apicaldata[[i]])
  anova_output_HR[i] <- summary(temp2)
  rm(temp2)
}
names(anova_output_HR) <- chemnames
print(anova_output_HR)
anova_HR_summary <- plyr::ldply(anova_output_HR, broom::tidy) %>%
  mutate(adj.p.value = p.adjust(p.value, "fdr"))
write_csv(anova_HR_summary, file = "anova_HR_summary.csv")

#Affected Rate
anova_output_AR <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- aov(Dose..mg.L. ~ Affected.Rate, data = apicaldata[[i]])
  anova_output_AR[i] <- summary(temp2)
  rm(temp2)
}
names(anova_output_AR) <- chemnames
print(anova_output_AR)
anova_AR_summary <- plyr::ldply(anova_output_AR, broom::tidy) %>%
  mutate(adj.p.value = p.adjust(p.value, "fdr"))
write_csv(anova_AR_summary, file = "anova_AR_summary.csv")

####Means and SD####
apicaldatasummary <- list()
for(i in 1:length(apicaldata)) {
  temp2 <- apicaldata[[i]] %>% group_by(Dose..mg.L.) %>%
    summarise("Survival Rate" = mean(Survival.Rate, na.rm=TRUE), SD_SR = sd(Survival.Rate, na.rm=TRUE),
              "Deformity Rate" = mean(Deformity.Rate, na.rm=TRUE), SD_DR = sd(Deformity.Rate, na.rm=TRUE),
              "Hatch Rate" = mean(Hatch.Rate, na.rm=TRUE), SD_HR = sd(Hatch.Rate, na.rm=TRUE),
              "Affected Rate" = mean(Affected.Rate, na.rm=TRUE), SD_AR = sd(Affected.Rate, na.rm=TRUE))
  colnames(temp2) <- c("Dose", "Survival_Rate", "SR_SD", "Deformity_Rate", "DR_SD", "Hatch_Rate", "HR_SD", "Affected_Rate", "AR_SD")
  apicaldatasummary[[i]] <- as.data.frame(temp2)
  rm(temp2)
}
names(apicaldatasummary) <- chemnames
print(apicaldatasummary)
apicaldatasummary_tibble <- apicaldatasummary %>%
  plyr::ldply()

#Output Sheets
all_apicaldata <- plyr::ldply(apicaldata) %>%
  dplyr::select(-.id)
write_csv(all_apicaldata, file = "master_apical_data.csv")

####Levene Test####
# Homogeneity of variance -- of Affected Rate
LeveneResults <- all_apicaldata %>%
  group_by(Chemical) %>%
  summarise(car::leveneTest(Affected.Rate, as.factor(Dose..mg.L.)))
#Adding an is.significant column to easily parse significant values in a spreadsheet
LeveneResults <- LeveneResults %>%
  na.omit() %>%
  mutate(is.significant = if_else(
    condition = `Pr(>F)` < 0.05,
    true = TRUE,
    false = FALSE
  ))
DT::datatable(LeveneResults)

write_csv(x = LeveneResults, file = "Levene_Test_Results.csv")
#rm(VarianceCheck)


####Dunnett's Test####
#PostHoc tests
set.seed(2345)

#Dunnett's Test
Dunnett_results <- apicaldata_ol_rm %>%
  group_by(Chemical) %>%
  nest() %>%
  mutate(model = map(data, ~ DunnettTest(
    x = .$Affected.Rate, g = .$Dose..mg.L.
  ), data = .)) #Performing the Dunnett's test and saving it is a variable

#Creating list of summaries
#Since the PostHocTest object cannot be coerced to a tidy tibble using broom::tidy()... we got creative
Dunnett_list <-
  list() #What we are trying to do is index the results and see what the significant results were... so we are using a list which can be later coerced into a tibble to easily index...
for (i in 1:length(unique(apicaldata_ol_rm$Chemical))) {
  Dunnett_list[[Dunnett_results$Chemical[i]]] <-
    Dunnett_results$model[[i]][["0"]] %>% #Take Dunnett's test results without any of the fancy summary information and shove it into a named list
    as.data.frame() %>% #Coerce to a data frame temporarily so what we can take the row names of the reults and turn them into a variable with rownames_to_column
    rownames_to_column(var = "dose")
}
Dunnett_comb <-
  plyr::ldply(Dunnett_list) #this function combines all of the lists together and gives them a variable name according to the chemical
Dunnett_comb$dose = substr(Dunnett_comb$dose,
                           start = 1,
                           stop = nchar(Dunnett_comb$dose) - 2) %>%
  as.numeric() #Here we are fixing the dose column... the dose column has the test dose related to the control... but we just want to see what the test dose is without it giving us redundant information about the comparison to the control for every observation...
Dunnett_comb <- as_tibble(Dunnett_comb) #Coerce to a tidy tibble
#Great, a nice tibble that we can export

#Now just to add one more column
Dunnett_comb <- Dunnett_comb %>%
  mutate(adj_p.value = p.adjust(pval, method = "fdr")) %>%
  mutate(is.significant = if_else(
    condition = adj_p.value < 0.05,
    true = TRUE,
    false = FALSE
  )) %>%
  rename(Dose = dose)
DT::datatable(Dunnett_comb)

#Write to a .csv
write_csv(Dunnett_comb, file = "Dunnett_test_results.csv")

# #indexing what the significant results were...
Dunnett_Sig_Results <-
  Dunnett_comb[which(Dunnett_comb$adj_p.value <= 0.05), ]
Dunnett_Sig_Results
# #Cool!

##################

#Survival Rate
#Saving Images
for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
  ggplot(aes(x=Dose, y=Survival_Rate)) +
    geom_point(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Survival_Rate-SR_SD, ymax=Survival_Rate+SR_SD), width=.2,
                  position=position_dodge(.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    scale_y_continuous(name = "Survival Rate",
                       breaks = pretty_breaks(),
                       expand = expansion(mult = (c(0, 0))),
                       limits = c(0,1.2))
                       
    print(p)
    ggsave(paste0(names(apicaldatasummary[i]), ".png"), device = "png", path = paste0(getwd(), "/Images/Survival_Rate"))
}

#Deformity Rate
for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
    ggplot(aes(x=Dose, y=Deformity_Rate)) +
    geom_point(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Deformity_Rate-DR_SD, ymax=Deformity_Rate+DR_SD), width=.2,
                  position=position_dodge(.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    scale_y_continuous(name = "Rate of Deformity",
                       breaks = pretty_breaks(),
                       expand = expansion(mult = (c(0, 0))),
                       limits = c(-0.1,1.2))
  print(p)
  ggsave(paste0(names(apicaldatasummary[i]), ".png"), device = "png", path = paste0(getwd(), "/Images/Deformity_Rate"))
}

#Hatch Rate
for(i in 1:length(apicaldatasummary)) {
  p <- apicaldatasummary[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
    ggplot(aes(x=Dose, y=Hatch_Rate)) +
    geom_point(stat="identity", color="black", fill="white", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Hatch_Rate-HR_SD, ymax=Hatch_Rate+HR_SD), width=.2,
                  position=position_dodge(.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title=names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    scale_y_continuous(name = "Hatch Rate",
                       breaks = pretty_breaks(),
                       expand = expansion(mult = (c(0, 0))),
                       limits = c(-0.1,1.5))
  print(p)
  ggsave(paste0(names(apicaldatasummary[i]), ".png"), device = "png", path = paste0(getwd(), "/Images/Hatch_Rate"))
}

#Affected Rate
#Saving Images
apicaldatasummary_2 <- inner_join(apicaldatasummary_tibble, Dunnett_comb, by = c("Dose", ".id")) %>%
  group_by(.id) %>%
  split(f = .$.id)

for(i in 1:length(apicaldatasummary_2)) {
  p <- apicaldatasummary_2[[i]] %>%
    as.data.frame() %>%
    mutate(Dose = as.factor(Dose)) %>%
    ggplot(aes(x = Dose, y = Affected_Rate)) +
    geom_point(
      stat = "identity",
      color = "black",
      fill = "white",
      position = position_dodge()
    ) +
    geom_errorbar(
      aes(ymin = Affected_Rate - AR_SD, ymax = Affected_Rate + AR_SD),
      width = .2,
      position = position_dodge(.9)
    ) +
    geom_text(
      aes(
        label = if_else(
          condition = adj_p.value > 0.1,
          true = "",
          false = if_else(
            condition = adj_p.value <= 0.1 &
              adj_p.value > 0.05,
            true = "",
            if_else(
              condition = adj_p.value <= 0.05 &
                adj_p.value > 0.01,
              true = "*",
              false = if_else(
                condition = adj_p.value <= 0.01 &
                  adj_p.value > 0.001,
                true = "**",
                false = if_else(
                  adj_p.value <= 0.001 &
                    adj_p.value >= 0,
                  true = "***",
                  false = ""
                )
              )
            )
          )
        ),
        group = Dose,
        y = Affected_Rate + AR_SD + 0.1,
        x = Dose
      ),
      color = "red",
      size = 3.75
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    labs(title = names(apicaldatasummary)[i]) +
    xlab("Dose (mg/L)") +
    ylab("Affected Rate") +
    scale_y_continuous(
      name = "Affected Rate",
      breaks = pretty_breaks(),
      expand = expansion(mult = (c(0, 0))),
      limits = c(-0.1, 1.2)
    )
  print(p)
  ggsave(
    paste0(names(apicaldatasummary[i]), ".png"),
    device = "png",
    path = paste0(getwd(), "/Images/Affected_Rate")
  )
}
