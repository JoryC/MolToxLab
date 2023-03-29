library(dplyr)

# You must compile BMR results before you can compile SS resutlss
#Please use runCompile_BMR_Results.sh to compile BMR results. It is a bash script that sources this R script in a loop

args <- commandArgs() #list the arguments that were used on command line to open the script
print(args)

endpoint_R <- c("negative", "positive")
direction <- c("neg", "pos")
step <- c("3_1-", "3_2-")

endpoint <- endpoint_R[as.numeric(args[6])] #index which endpoint to use for this iteration of the loop after running the script with `runCompile_BMR_Results.sh`
dir <- direction[as.numeric(args[6])]
step <- step[as.numeric(args[6])]
k=endpoint


sample_size <- 1000 #Change this to whatever sample size you want. Test with 10 or 100, get results with 1000 samples
n_dose_groups <- 6 #Change this to the number of dose groups that you have (including the control group)


rmarkdown::render(
  input = paste0(step, "CurveP_estimate_bmr_z-score_", k, "_direction", ".Rmd"),
  output_file = paste0("report_", step, "CurveP_estimate_bmr_z-score_", k, "_direction_", Sys.Date(), ".html"),
  rmarkdown::html_document(
    toc = TRUE,
    toc_depth = 5,
    toc_float = TRUE,
    number_sections = TRUE,
    code_folding = "show",
    df_print = "paged",
    code_download = TRUE,
    theme = "readable"
  ),
  params = list(
    direction = dir,
    sample_size = sample_size,
    n_dose_groups = n_dose_groups
  )
)
print(paste0("DONE!!!!!_", "report_", step, "CurveP_estimate_bmr_z-score_", k, "_direction_", Sys.Date(), ".html"))
#Free up memory
rm(list = ls())
gc(reset = TRUE, full = TRUE)
#rstudioapi::executeCommand("restartR")
q(save = "no")