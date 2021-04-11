#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to calculate the admixture of the individuals from the filtered and thinned hard-called data.
## to run the admixture software we first need to create .bed files from the .vcf file. We
## also need to change chromosome names to numbers in the .bim files from the .vcf file.

module load plink # Load module needed.

plink --vcf <Path_and_File_Name_of_VCF_File> --allow-extra-chr --make-bed --out <Output_Path_and_File_Name_1> # Create the .bed file needed.

awk '{$1=0 ; print ;}' <Path_and_File_Name_of_BIM_File> > <Output_Path_and_File_Name_1> # Change the chromosome name to number.

for K in `seq 1 20`;
    do admixture --cv <Path_and_File_Name_of_BED_File> $K | tee <Output_Path>/log${K}.out; # Calculate admixture for K 1 to 20.
done