#!/bin/bash

## Carla Coll Costa, C3.

## Script to calculate the depth or coverage of the data. 
## The script reads a .txt file in which each row contains information on one sample. The first column
## has the sample name and the second column the path to the .bam file of that sample. 

module load biokit # Load module needed.

while read sample path; do
    samtools depth $path |  awk '{sum+=$3} END {print sum/NR}' >> <Output_Path_and_File_Name>; # Calculate the depth of the sample and output it.
done < <Input_File_to_Read> # Give the file to read line by line.

paste <Input_File_to_Read> <Output_Path_and_File_Name> >> <Name_Temporary_File.tab> # Create a temporary tab file with the information of the input file to read ant the depth of that sample.
awk '{print $1"\t"$3}' <Name_Temporary_File.tab> >> <Final_File_with_Depths.tab> # Create the final file with each row containing the name of the sample and the mean depth.
rm <Output_Path_and_File_Name> # Remove temporary files created during the analysis.
rm <Name_Temporary_File.tab> # Remove temporary files created during the analysis.