Order of pipeline:
1. Knit Import_and_Exploratory_Data_Analysis from .rmd
	1.1. Images will be in ./Output/Images and speadsheets in ./Output
2. Knit transform raw data to z_scores from .rmd
	2.1. Image will be in ./Output/Images and spreadsheets in ./Output
3. Execute runCompile_BMR_Results.sh in bash terminal
4. Execute runCompile_SS_Results.sh in bash terminal
5. Results are in Output folder and Outpur/BMDs folder
