#### Libraries ####
library(Rcurvep)
library(tidyverse)
library(future)
library(future.apply)
library(furrr)
library(kableExtra)
source("Functions/RCurveP_Data_Wrangling_Functions.R")

#### calculating BMR ####

#read z-score dat
z_dat <- read_csv(file = "Data/z-score_data.csv")
#test <- readRDS(file = "~/MolToxLab/Behaviour/Jory_Behaviour_R_Markdown/Data/simi_norm_tib_4rcurvep_pearson.rds")

#metadata
metaData <- read.csv(file = "Data/MetaData.csv") #Import the Meta Data that includes information about the data in the folders
#CAS is the Chemical Abstract Service, MOA is the Mode of Action. This table includes useful information about the exposure concentrations for each chemical dose in mg/L. We'll use this later to create our final data frame
doseData <- metaData %>%
  select(plate_id, Dose1:Control) %>%
  gather(key = Dose, value = "Dose_mg_L", Dose1:Control)
number_of_dose_groups_per_chem <- 6 # Set this to the number of dose groups you have (including control)
lowest_dose_treatment_group <- doseData %>%
  filter(Dose != "Control") %>%
  dplyr::summarize(lowest_dose_treatment_group = min(Dose_mg_L))
lowest_dose_treatment_group <- lowest_dose_treatment_group[1,][[1]]
dir <- "pos " #Direction

z_dat <- z_dat %>%
  mutate(endpoint = "z_score") %>%
  rename(resp = z_score,
         chemical = Chemical) %>%
  mutate(`Dose(mg/L)` = `Dose(mg/L)`/lowest_dose_treatment_group[[1]]) %>%
  mutate(conc = log10(`Dose(mg/L)`)) %>%
  group_by(chemical, endpoint) %>%
  mutate(control_conc = max(conc)-5) %>% #making sure there are no infinite values for conncentrations. Controls are made to be 1 dose group down from the lowest concentration
  mutate(conc = if_else(condition = conc == -Inf, true = control_conc, false = conc)) %>%
  select(endpoint, chemical, conc, resp)

set.seed(550)

#sample size
ss <- 10

#BMR estimating step 1 in the negative direction
BMR_training_set <- expression(
  value(
    future({
      combi_run_rcurvep(
        z_dat,
        n_samples = ss, #Increase this number to 100 or 1000 for better results (takes a long time to run)
        keep_sets = c("act_set", "resp_set", "fp_set"),
        TRSH = seq(0.25, 6, by = 0.25),
        RNGE = 1000000,
        TrustHi = FALSE #We assume that there is more noise at the higher concentration, so if two correction sets are equal, go with the lower concentration set
      )
    }, seed = 550)
  )
)

BMR_training_set_act <- eval(BMR_training_set)

#Bmr estimation step 2
bmr_output <-
  estimate_dataset_bmr(BMR_training_set_act, plot = FALSE)
bmr_output_outcome <- bmr_output$outcome

plot(bmr_output)

#Outcome of BMR estimation
bmr_output_outcome %>%
  select(endpoint, qc, cor_exp_fit, cor_lm_fit, bmr_exp, bmr_ori) %>%
  kable(
    col.names = c(
      "Endpoint",
      "Quality Control Message",
      "Correlation of expotential fit",
      "Correlation of linear model fit",
      "BMR of exponential model",
      "BMR of linear model"
    ),
    align = 'llrrrr',
    caption = "Table 8. Summary of the estimated BMR"
  ) %>%
  kable_styling() %>%
  scroll_box()

#Pull out the BMR
bmr_output$outcome[which(is.na(bmr_output$outcome))] = 0
bmr_thresh <- bmr_output$outcome %>%
  mutate(bmr = if_else(condition = .$cor_lm_fit > .$cor_exp_fit, true = .$bmr_ori, false = .$bmr_exp)) %>%
  pull(bmr)

#Create input tibble for Benchmark Dose Calculation
input_tib <- bmr_output_outcome %>%
  nest_join(z_dat, by = c("endpoint"), keep = TRUE, name = "data") %>%
  select(RNGE, endpoint, bmr_exp, data)

#calculate the BMD step 1
system.time(
  bmd_dat_neg <-
    future_pmap(
      input_tib,
      ~ combi_run_rcurvep(
        ..4,
        TRSH = ..3,
        RNGE = ..1,
        n_samples = ss,
        keep_sets = c("act_set", "resp_set"),
        TrustHi = FALSE
      ),
      .options = furrr_options(seed = 550)
    )
)

#calculate the BMD step 2
sum_bmd_dat_neg <-
  summarize_rcurvep_output(
    bmd_dat_neg[[1]],
    inactivate = "INVERSE",
    ci_level = 0.95
  )

#summarise the results
simulated_curves_summary <-
  summarise_rcurvep_results(
    resp_set = sum_bmd_dat_neg$result$resp_set,
    act_set = sum_bmd_dat_neg$result$act_set,
    act_summary = sum_bmd_dat_neg$act_summary,
    reject_hit_conf_under = 0.5
  )

bmds <- bmd_results(simulated_curves_summary)

curve_hits <- simulated_curves_summary %>%
  filter(hit == 1)

confident_curves <- confident_hits(summary_dat = simulated_curves_summary, reject_hit_conf_under = 0.5)
# If Error says "out of bounds", reduce conf filter

p <- ggplot(NULL, aes(conc, resp)) +
  geom_point(
    data = simulated_curves_summary,
    mapping = aes(x = conc, y  = resp, group = sample_id),
    color = "gray75",
    alpha = 0.2
  ) +
  geom_point(
    data = curve_hits,
    mapping = aes(x = conc, y = resp, group = sample_id),
    color = "#a9a9a9",
    alpha = 0.7
  ) +
  geom_line(
    data = simulated_curves_summary,
    mapping = aes(x = conc, y = resp, group = sample_id),
    color = "gray75",
    alpha = 0.2
  ) +
  geom_line(
    data = curve_hits,
    mapping = aes(x = conc, y = resp, group = sample_id),
    color = "#a9a9a9",
    alpha = 0.7
  ) +
  # geom_line(data = confident_curves[["cil"]], mapping = aes(x = conc, y = resp, group = sample_id), color = "purple") +
  # geom_line(data = confident_curves[["ciu"]], mapping = aes(x = conc, y = resp, group = sample_id), color = "yellow") +
  # geom_line(data = confident_curves[["mean"]], mapping = aes(x = conc, y = resp, group = sample_id), color = "green") +
  # geom_line(data = confident_curves[["median"]], mapping = aes(x = conc, y = resp, group = sample_id), color = "blue") +
  geom_hline(
    yintercept = bmr_thresh,
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(data = bmds,
             aes(xintercept = median_POD),
             color = "blue") +
  geom_vline(data = bmds, aes(xintercept = ciu_POD), color = "yellow") +
  # geom_vline(data = bmds_pearson, aes(xintercept = mean_POD), color = "green") +
  geom_vline(data = bmds, aes(xintercept = cil_POD), color = "purple") +
  geom_text(
    data = bmds,
    mapping = aes(label = paste("Hit =", hit_confidence), x = lowest_conc +
                    1),
    y = -85,
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap( ~ chemical, scales = "free") +
  ylim(-6, 6) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey50", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#save the image
ggsave(
  plot = p,
  filename = paste0(
    "Output/Images/", "z_score_", dir, "_", "_processed_curves_",
    ss,
    "_samples.svg"
  ),
  device = "svg",
  units = "px",
  width = 1920,
  height = 1080,
  bg = "white",
  scale = 2
)

#Distribution of BMDs
p_distribution_bmds <- curve_hits %>%
  group_by(chemical) %>%
  ggplot(data = ., aes(x = POD)) +
  geom_histogram() +
  facet_wrap( ~ chemical, scales = "free_x") +
  geom_vline(aes(xintercept = median_POD), color = "blue") +
  geom_vline(aes(xintercept = cil_POD), color = "purple") +
  geom_vline(aes(xintercept = ciu_POD), color = "yellow") +
  geom_vline(aes(xintercept = lowest_conc),
             color = "black",
             linetype = "dashed") +
  geom_vline(aes(xintercept = highest_conc),
             color = "black",
             linetype = "solid") +
  theme_classic() +
  scale_x_continuous(expand = expansion(add = 1)) +
  xlab("BMD") +
  geom_text(
    mapping = aes(label = paste("Hit =", hit_confidence), x = lowest_conc +
                    2),
    y = 1000,
    inherit.aes = FALSE,
    size = 4
  )

#save the image
ggsave(
  plot = p_distribution_bmds,
  filename = paste0(
    "Output/Images/", "z_score_", dir, "_", "_distribution_bmds_",
    ss,
    "_samples.svg"
  ),
  device = "svg",
  units = "px",
  width = 1920,
  height = 1080,
  bg = "white",
  scale = 2
)

#Table the results
temp <- list(bmds)
names(temp) <- c(paste0("z_score_", "positive_dir"))
all_summarized_bmd_dat <- plyr::ldply(temp)
all_summarized_bmd_dat <- all_summarized_bmd_dat %>%
  arrange(chemical)

all_summarized_bmd_dat %>%
  mutate(lowest_conc = lowest_conc + 1) %>%
  mutate(
    lowest_conc = (10 ^ lowest_conc)*lowest_dose_treatment_group,
    highest_conc = sprintf((10 ^ highest_conc)*lowest_dose_treatment_group, fmt = '%.7f'),
    cil_POD = (10 ^ cil_POD)*lowest_dose_treatment_group,
    POD_med = (10 ^ POD_med)*lowest_dose_treatment_group,
    ciu_POD = (10 ^ ciu_POD)*lowest_dose_treatment_group,
    median_POD = (10 ^ median_POD)*lowest_dose_treatment_group,
    mean_POD = (10 ^ mean_POD)*lowest_dose_treatment_group
  ) %>%
  mutate(POD_med = if_else(is.na(POD_med), true = median_POD, false = POD_med)) %>%
  select(
    chemical,
    .id,
    lowest_conc,
    highest_conc,
    cil_POD,
    POD_med,
    ciu_POD,
    hit_confidence
  ) %>%
  kable(
    col.names = c(
      "Chemical",
      "endpoint",
      "Low Dose (mg/L)",
      "High Dose (mg/L)",
      "BMD - low confidence interval (0.05)",
      "BMD",
      "BMD - high confidence interval (0.95)",
      "Hit Confidence"
    ),
    caption = "Table 9. Summary of all the benchmark doses (BMDs) of every chemical"
  ) %>%
  kable_styling() %>%
  row_spec(
    which(
      all_summarized_bmd_dat$median_POD < all_summarized_bmd_dat$highest_conc
    ),
    bold = TRUE,
    background = "gray"
  ) %>%
  scroll_box()

#Save the data
data <- all_summarized_bmd_dat %>%
  mutate(lowest_conc = lowest_conc + 1) %>%
  mutate(
    lowest_conc = (10 ^ lowest_conc)*lowest_dose_treatment_group,
    highest_conc = sprintf((10 ^ highest_conc)*lowest_dose_treatment_group, fmt = '%.7f'),
    cil_POD = (10 ^ cil_POD)*lowest_dose_treatment_group,
    POD_med = (10 ^ POD_med)*lowest_dose_treatment_group,
    ciu_POD = (10 ^ ciu_POD)*lowest_dose_treatment_group,
    median_POD = (10 ^ median_POD)*lowest_dose_treatment_group,
    mean_POD = (10 ^ mean_POD)*lowest_dose_treatment_group
  ) %>%
  mutate(POD_med = if_else(is.na(POD_med), true = median_POD, false = POD_med))

write_csv(x = data, file = paste0("Output/", "BMDs/", "z_score_", dir, "_", "BMDs_", ss, "_samples.csv"))
