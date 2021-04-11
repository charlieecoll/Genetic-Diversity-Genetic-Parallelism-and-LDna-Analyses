#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to retrieve the 2 dimensional site frequency sprectrum of each pair of populations using genotype 
## likelihoods. The script is automated in arrays so each array will do the calculations for one 
## pair of populations. Each array will read a line of a file. Each line contains information on one of the pairwise 
## comparisons: each of the two columns has the path to the .saf.idx file generated for one of the populations
## during the calculations of the 1 dimensional site frequency spectrum. There are as many rows as pairwise
## comparisons.

conda activate angsd0921 # Load module needed.

sed -n ${SLURM_ARRAY_TASK_ID}p <Path_and_File_to_be_Read_by_Arrays> > tmp_${SLURM_ARRAY_TASK_ID} # Set temporary files for the array.

while read pop1 pop2; do
    first=`echo ${pop1} | sed 's|.*/||' | sed 's/\..*//'` # Get the name of the population in the first column.
    second=`echo ${pop2} | sed 's|.*/||' | sed 's/\..*//'` # Get the name of the population in the second column.
    realSFS ${pop1} ${pop2} > <Output_Path>/${first}_${second}.sfs # Get the 2D SFS.
done < tmp_${SLURM_ARRAY_TASK_ID} # Each array reads one line of the input file.

rm tmp_${SLURM_ARRAY_TASK_ID} # Remove temporary files created during the analysis.
