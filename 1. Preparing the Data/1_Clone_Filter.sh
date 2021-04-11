#!/bin/bash

## Carla Coll Costa, C3.

## Script to remove the PCR duplicates with clone_filter.
## The script reads a .txt file with several rows and columns. Each row contains data on one sample 
## (or fq.gz file). The first column corresponds to a unique ID of the sample (samples might have
## been sequenced more than once in which case we will have more than one file for the same sample).
## The second column is the sample name (same name for samples that have been resequenced), the third
## column has the population name of the sample and the fourth the path to the fq.gz files.

module load stacks # Load module needed.

while read -r unique_id sample population path; do
    mkdir -p /${population}/${sample}/ # Create a new directory for each of the populations and each of the samples.
    output_directory="/${population}/${sample}/" # Choose the directory created as the output directory.
    clone_filter -1 <Path_to_fq.gz_File>/${population}/${unique_id}_1.fq.gz -2 <Path_to_fq.gz_File>/${population}/${unique_id}_2.fq.gz -i gzfastq -o $output_directory # Apply clone filter giving as input the 2 fq.gz files per sample (forward and reverse) and output the files in the output directory chosen.
done <  <Input_File_to_Read> # Give the file to read line by line.
