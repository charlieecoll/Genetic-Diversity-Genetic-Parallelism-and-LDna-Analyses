#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to calculate the runs of homozygosity from the filtered hard called data.
## To calculate ROHs it is necessary to have the allele frequencies in the vcf file that is going to be used
## and we can create a file with allele frequencies from that vcf file. At the end it is necessary to create
## a map file from the ROH file to be able to read the data in RStudio.

module load samtools # Load modules needed.

bcftools +fill-tags <Path_and_File_Name_of_VCF_File> > <Output_Path_and_Name_of_VCF_File_with_AF> # Add all possible tags (including AC and AF) to the vcf file.

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%AF\n' <Output_Path_and_Name_of_VCF_File_with_AF> | bgzip -c > <Output_Path_and_File_Name_of_AF_File.tab.gz> && tabix -s1 -b2 -e2 <Output_Path_and_File_Name_of_AF_File.tab.gz> # Create a file with allele frequencies from the VCF file with AF.

bcftools roh --AF-file <Output_Path_and_File_Name_of_AF_File.tab.gz> --samples-file <TXT_File_with_Samples_to_Include_(samplenameperline)> --output <Output_Path_and_File_Name_with_ROH_Data> <Input_VCF_File_with_AF_tag>

module load plink # Load modules needed.

plink --vcf <Input_VCF_File_with_AF_tag> --allow-extra-chr --recode --out <Output_Path_and_File_Name_of_MAP_File> # Create mapfile needed to process lots of ROH data in Rstudio