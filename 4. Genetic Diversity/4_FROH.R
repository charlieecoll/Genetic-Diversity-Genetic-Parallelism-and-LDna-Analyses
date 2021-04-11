#################################################################
##                             ROH                             ##
#################################################################

## Carla Coll Costa, C3. January 2021.
## This script takes the files on the ROH data calculated in BCFTOOLS with the script 
## 2_ROH. It also needs a map file and a genotype file for the ROH data calculated. 


# Upload the required packaged
library(detectRUNS)
library(stringr)
library(Eagle)
library(easyGgplot2)
library(plyr)

# Set the working directory
setwd("<Working_Directory>")

# Take the data
data = readExternalRuns(inputFile="<Input_File_with_ROH_Data>", program="BCFtools")

# Remove the scaffolds
data = data[- grep("scaffold", data$chrom),]

# Get the map path file
map = file.path("<Input_Map_File>")

# Get the genotype file
genotype=file.path("~/Desktop/Master's Thesis/Results/17. ROH/0_data/Genotype_File.ped")

# We need to add the population to the group column
data$group = data$id
data$group = sub("-[^-]+$", "", data$group)

# We substitute the chromosome name to an appropriate name
data$chrom = str_replace_all(data$chrom,"groupI\\b","1")
data$chrom = str_replace_all(data$chrom,"groupII\\b","2")
data$chrom = str_replace_all(data$chrom,"groupIII\\b","3")
data$chrom = str_replace_all(data$chrom,"groupIV\\b","4")
data$chrom = str_replace_all(data$chrom,"groupV\\b","5")
data$chrom = str_replace_all(data$chrom,"groupVI\\b","6")
data$chrom = str_replace_all(data$chrom,"groupVII\\b","7")
data$chrom = str_replace_all(data$chrom,"groupVIII\\b","8")
data$chrom = str_replace_all(data$chrom,"groupIX\\b","9")
data$chrom = str_replace_all(data$chrom,"groupX\\b","10")
data$chrom = str_replace_all(data$chrom,"groupXI\\b","11")
data$chrom = str_replace_all(data$chrom,"groupXII\\b","12")
data$chrom = str_replace_all(data$chrom,"groupXIII\\b","13")
data$chrom = str_replace_all(data$chrom,"groupXIV\\b","14")
data$chrom = str_replace_all(data$chrom,"groupXV\\b","15")
data$chrom = str_replace_all(data$chrom,"groupXVI\\b","16")
data$chrom = str_replace_all(data$chrom,"groupXVII\\b","17")
data$chrom = str_replace_all(data$chrom,"groupXVIII\\b","18")
data$chrom = str_replace_all(data$chrom,"groupXIX\\b","19")
data$chrom = str_replace_all(data$chrom,"groupXX\\b","20")
data$chrom = str_replace_all(data$chrom,"groupXXI\\b","21")

# Get information on the populations and individuals
pop_labs = read.table("<File_with_Information_on_Population_and_Individuals>", header=TRUE)

#### FROH: ROH based inbreeding ####

# Calculate inbreeding
inbreeding = Froh_inbreeding(runs=data, mapFile=map, genome_wide=TRUE)

# Merge other information to plot latter
inbreeding = merge(inbreeding, pop_labs, by.x="id", by.y="Individual")

# Delete duplicated column
inbreeding$group = NULL

# Output a tab separated file with th information on inbreeding based on runs of homozygosity
write.table(inbreeding, "<Output_File_Name>", sep="\t")
