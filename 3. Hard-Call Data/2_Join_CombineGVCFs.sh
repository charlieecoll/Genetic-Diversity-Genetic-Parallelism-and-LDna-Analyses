#!/bin/bash

## Carla Coll Costa, C3.

## Script to combine the gvcf files for each of the samples in one gvcf file.
## The script expects the user to give the path and name of the GVCF file for each of the samples
## just before the --variant option. (As many "--variant" as files there are).

module load gatk # Load module needed.

gatk CombineGVCFs \
-R <Reference_Genome_with_Path> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
--variant <Path_and_Name_of_GVCF_File> \
-O <Output_Path_and_File_Name>
