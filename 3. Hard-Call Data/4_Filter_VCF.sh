#!/bin/bash

## Carla Coll Costa, C3.

## Script to filter and thin the called data according to different filtering options.

module load vcftools # Load module needed.

vcftools --vcf <Path_and_Name_of_Input_VCF_File> --keep IndDepth2 --max-alleles 2 --min-alleles 2 --minGQ 20 --max-missing 0.9 --maf 0.01 --mac 2 --remove-indels --recode --out <Path_and_Output_File_Name.vcf> # Filter the data.

vcftools --vcf <Path_and_Name_of_Input_VCF_File> --thin 10000 --recode --out <Path_and_Output_File_Name> # Thin the data.
