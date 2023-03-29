library(dplyr)

#Before Compiling_SS_Results, you first must compile BMR results using runCompile_BMR_results.sh
#Please use runCompile_SS_Results.sh to compile all of the results. It is a bash script that sources this R script in a loop

args <- commandArgs() #list the arguments that were used on command line to open the script
print(args)

endpoint_R <- c("negative", "positive")
method_R <- c("CurveP", "Fitted")
direction <- c("neg", "pos")
step <- c("4_1-", "4_2-")

endpoint <- endpoint_R[as.numeric(args[6])] #index which endpoint to use for this iteration of the loop after running the script with `runCompile_BMR_Results.sh`
dir <- direction[as.numeric(args[6])]
step <- step[as.numeric(args[6])]
method <- method_R[as.numeric(args[7])] #index which method to use for this iteration of the loop after running the script with `runCompile_SS_Results.sh`
k=endpoint
i=method

n_dose_groups <- 6 #Change this to the number of dose groups that you have (including the control group)
sample_size <- 1000 #Change this to whatever sample size you want. Test with 10 or 100, get results with 1000 samples

rmarkdown::render(
  input = paste0(step, i, "_bmd_z-score_", k, "_direction", ".Rmd"),
  output_file = paste0("report_", step, i, "_bmd_z-score_", k, "_direction", ".html"),
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
print(paste0("DONE!!!!!_", "report_", step, i, "_bmd_z-score_", k, "_direction", ".html"))
#Free up memory
rm(list = ls())
gc(reset = TRUE, full = TRUE)
#rstudioapi::executeCommand("restartR")
q(save = "no")
