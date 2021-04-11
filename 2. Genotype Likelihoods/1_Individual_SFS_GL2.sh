#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to get the SNP variants of each individual using genotype likelihoods. 
## This script reads a file where each line corresponds to an individual. The first column has the name
## of the individual and the second column the path to the .bam file of that individual.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -C 50" # Filters applied in the analysis.
    # GENERAL FOR ANGSD.
    # -uniqueOnly 1 : Remove reads that have multiple best hits.
    # -remove_bads 1: Removes reads with a flag above 255, this is, the secondary alignments of some reads.
    # -minMapQ 20: Minimum MapQ quality.
    # -minQ 20: Minimum base quality score.
    # -C 50: Adjust mapQ for excessive mismatches (50).
TODO="-dosaf 1 -GL 2" # Filters applied in the analysis.
    # CONCRETE FOR SFS.
    # -dosaf 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE.
    # -GL 1: Genotype likelihood model.

# -P: Number of threads. 
# -fold 1: If you don't have the ancestral state, you can instead estimate the folded SFS. This is done by supplying the -anc with the reference genome and applying -fold 1 to realSFS.

conda activate angsd0921 # Load module needed.

# ANGSD: SFS BASED ON GENOTYPE LIKELIHOODS
while read sample_name path; do 
    angsd  $FILTERS $TODO -ref <Reference_Genome_with_Path> -anc <Reference_Genome_with_Path> -i $path -P 16  -out <Output_Path>/${sample_name};
done < <Input_File_to_Read> # Give the file to read line by line.
# Output: 1) angsdput.saf.idx 2) angsdput.saf.pos.gz 3) angsdput.saf.gz

# realSFS: Get ML estimates. This program will estimate the (multi) SFS based on a .saf file generated from the ./angsd [options] -doSaf .
while read sample_name path; do 
    realSFS <Output_Path>/${sample_name}.saf.idx > <Output_Path>/${sample_name}.sfs;
done < <Input_File_to_Read> # Give the file to read line by line.

