library(dplyr)

# You must compile BMR results before you can compile SS resutlss
#Please use runCompile_BMR_Results.sh to compile BMR results. It is a bash script that sources this R script in a loop

args <- commandArgs() #list the arguments that were used on command line to open the script
print(args)

endpoint_R <- c("pearson", "spearman")

endpoint <- endpoint_R[as.numeric(args[6])] #index which endpoint to use for this iteration of the loop after running the script with `runCompile_BMR_Results.sh`
k=endpoint

value <- c("totaldist") # Can add other endpoints to this like lardist, smldist, etc.
sample_size <- 1000 #Change this to whatever sample size you want. Test with 10 or 100, get results with 1000 samples
method <- c("rcurvep")

rmarkdown::render(
  input = paste0("estimate_bmr_", k, "_rcurvep.rmd"),
  output_file = paste0("report_", "estimate_bmr_", k, "_rcurvep", "_", Sys.Date(), ".html"),
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
    endpoint = k,
    method = method
  )
)
print(paste0("Done", "report_", "estimate_bmr_", k, "_rcurvep", "_", Sys.Date(), ".html"))
#Free up memory
rm(list = ls())
gc(reset = TRUE, full = TRUE)
#rstudioapi::executeCommand("restartR")
q(save = "no")