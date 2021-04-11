#!/bin/bash

## Carla Coll Costa, C3.

## Script to call the genotypes from the joint GVCF file with information on all SNPs of all
## the samples (created in 10_Join_CombineGVCFs.sh).

module load gatk # Load module needed.

gatk GenotypeGVCFs -R <Reference_Genome_with_Path> -V <Join_GVCF_File_Created_with_CombineGVCFs> -O <Path_and_Output_File_Name.vcf>
