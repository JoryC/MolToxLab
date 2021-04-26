#### Generic script for cleaning .txt data from the Fluorometer in Mennigen Lab ####

#Plate Layout
	#Rows A-F all contain 9 fish (Row G is empty)
	#Row H has 4 fish in first 4 wells
#Use GNU Stream editor sed
sed -n -s '20d;14,21p' CHEMICAL_*.txt > CHEMICAL_Combined.csv
	#Delete line 20 (row G) and print out lines 14-21 from input (all replicate runs)  and export to output file (A merged .csv withh all/both replicates)
#Use cut and bash variable
badcol='1,2,12,13,14'
cut -f$badcol --complement CHEMICAL_Combined.csv
