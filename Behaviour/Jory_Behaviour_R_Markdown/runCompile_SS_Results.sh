#!/bin/bash

#Note: use runCompile_BMR_results.sh before runCompile_SS_results
#To ensure that Rstudio is being reset after each iteration of the loop to clear the memory, we will source the Compile R script using a loop in bash... Running .rs.Restart in an R for loop ends the loop, so this is my solution to that problem

for i in {1..2}
do
  for j in {1..2}
  do
  echo "Rscript Compile_SS_Results.R $i $j"
  Rscript Compile_SS_Results.R $i $j
  done
done