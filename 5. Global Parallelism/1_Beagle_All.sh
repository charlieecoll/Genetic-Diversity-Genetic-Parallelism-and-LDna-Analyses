#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to create a beagle file from all the .bam files with the information of the variants.
## The script expects a file in which each row has the path and file name to a .bam file of the 
## data set.

# GENERAL OPTIONS
# -GL 2: Genotype likelihood model. 
    # 2: Uses GATK.
# -nThreads: 
# -doGlf: Output the log genotype likelihoods to a file. 
    # 2: Beagle genotype likelihood format (use directly for imputation)
# -minInd: Only keep sites with at least minIndDepth (default is 1) from at least [int] individuals. 
    # 129: Half of the dataset.
# -skipTriallelic: Skip triallelic SNPs.
# -uniqueOnly 1: Remove reads that have multiple best hits.
# -remove_bads 1: Removes reads with a flag above 255, this is, the secondary alignments of some reads.
# -minMapQ 20: Minimum MapQ quality.
# -minQ 20: Minimum base quality score.
# -doMajorMinor: Many method assume that polymorphic sites are diallelic. For these methods one needs to define what is the major and minor allele. We allow the major and minor to be determined from either the counts of nucleotides, based on genotype likelihoods, specified by the ancestral/reference or even force both major minor to specific bases, which can be useful if you compare with HapMap data etc. 
    # 1: From input for either sequencing data like bam files or from genotype likelihood data like glfv3 the major and minor allele can be inferred directly from likelihoods. We use a maximum likelihood approach to choose the major and minor alleles.
# -doMaf: The major and minor allele is first inferred from the data. 
    # 2: Known major, Unknown minor. Here the major allele is assumed to be known (inferred or given by user) however the minor allele is not determined. Instead we sum over the 3 possible minor alleles weighted by their probabilities.
# -SNP_val: Only work with sites with a p-value less than [float].

conda activate angsd0921 # Load module needed.

angsd -GL 2 -out <Output_Path_and_File_Name> -nThreads 2 -doGlf 2 -minInd 129 -skipTriallelic -uniqueonly 1 -minMapQ 20 -minQ 20 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -bam <Input_File_with_Path_and_Name_to_BAM_Files> # Create one beagle file for all the .bam files.