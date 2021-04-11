#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to transform the vcf files in appropriate 012 files to calculate and plot the LD heatmap and do
## the LDna analysis.
## This script expects the vcf file with all the variants from all chromosomes and will create a 
## separate 012 output file for the chromosomes indicated from the vcf file. It also transforms the 
## -1 values of the 012 file to empty values.

vcftools --vcf <Input_Path_and_Name_of_VCF_File> --chr <Chromosome_to_Analyse> --012 --out <Output_Path_and_File_Name.012> # Keep the data from each of the chromosomes separately in a 012 file

sed -i "s,-1, NA, g" <Path_and_Name_of_012_File> # Transform de -1 of the .012 files in an empty value
