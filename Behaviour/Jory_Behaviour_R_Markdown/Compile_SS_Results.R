library(dplyr)

#Before Compiling_SS_Results, you first must compile BMR results using runCompile_BMR_results.sh
#Please use runCompile_SS_Results.sh to compile all of the results. It is a bash script that sources this R script in a loop

args <- commandArgs() #list the arguments that were used on command line to open the script
print(args)

endpoint_R <- c("pearson", "spearman")
method_R <- c("rcurvep", "fitted")

endpoint <- endpoint_R[as.numeric(args[6])] #index which endpoint to use for this iteration of the loop after running the script with `runCompile_SS_Results.sh`
method <- method_R[as.numeric(args[7])] #index which method to use for this iteration of the loop after running the script with `runCompile_SS_Results.sh`
i=endpoint
j=method

value <- c("totaldist") # Can add more endpoints to this like lardist, smldist, etc.
sample_size <- 100 #Change this to whatever sample size you want. Test with 10 or 100, get results with 1000 samples

rmarkdown::render(
  input = paste0(i, "_", j, ".rmd"),
  output_file = paste0("report_", i, "_", j, "_", Sys.Date(), ".html"),
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
    value = value,
    sample_size = sample_size,
    endpoint = i,
    method = j
  )
)
print(paste0("Done", "report_", i, "_", j, "_", Sys.Date(), ".html"))
#Free up memory
rm(list = ls())
gc(reset = TRUE, full = TRUE)
#rstudioapi::executeCommand("restartR")
q(save = "no")