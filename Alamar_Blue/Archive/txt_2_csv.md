#### Cleaning .txt data from the Fluorometer in Mennigen Lab ####

#Plate Layout assumption
	#Rows A-F all contain 9 fish (Row G is empty) [solved]
	#Row H has 4 fish in first 4 wells [define as NA in Rstudio.. work in progress]

#### Command Line/Terminal ####

#Delete line 20 (row G) and print out lines 14-21 from input (all replicate runs)  and export to output file (A merged .csv with all/both replicates)

#Use GNU Stream editor sed
	#For Baseline
sed -n -s '20d;14,21p' CHEMICAL_Baseline*.txt > CHEMICAL_Baseline_Combined.csv
	#For 24h
sed -n -s '20d;14,21p' CHEMICAL_24*.txt > CHEMICAL_24h_Combined.csv

#Define columns that are empty (1 & 2 residual from .txt file) and columns that contain no fish at all (columns 10, 11, and 12 on 96-well plate)
badcol='1,2,12,13,14' 

#Use cut and bash variable to cut and output one file input at a time
cut path/to/CHEMICAL_*_Combined.csv -f$badcol --complement > CHEMICAL_Cut.csv

#### In RStudio ####
