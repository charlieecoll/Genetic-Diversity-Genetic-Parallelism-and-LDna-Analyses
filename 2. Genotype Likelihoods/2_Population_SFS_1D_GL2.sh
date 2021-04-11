#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to retrieve the 1 dimensional site frequency sprectrum of each population using genotype 
## likelihoods. The script is automated in arrays so each array will do the calculations for one 
## population. Each array will read a line of a file with information on each population.
## The first column corresponds to the population name and the second column has the path to a file that contains
## the paths to the .bam files of the individuals of the corresponding population (each line has the path to 
## one .bam file).


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
    # -GL 2: Genotype likelihood model.

# -P: Number of threads. 
# -fold 1: If you don't have the ancestral state, you can instead estimate the folded SFS. This is done by supplying the -anc with the reference genome and applying -fold 1 to realSFS.

conda activate angsd0921 # Load module needed.

sed -n ${SLURM_ARRAY_TASK_ID}p <Path_and_File_to_Read_by_Each_Array> > tmp_${SLURM_ARRAY_TASK_ID} # Set temporary files for the array.

# ANGSD: SFS BASED ON GENOTYPE LIKELIHOODS
while read population bamfiles; do
    angsd  $FILTERS $TODO -ref <Reference_Genome_with_Path> -anc <Reference_Genome_with_Path> -bam $bamfiles  -out <Output_Path>/${population};
done < tmp_${SLURM_ARRAY_TASK_ID} # Each array reads one line of the input file.
# Output: 1) angsdput.saf.idx 2) angsdput.saf.pos.gz 3) angsdput.saf.gz

# realSFS: Get ML estimates. This program will estimate the (multi) SFS based on a .saf file generated from the ./angsd [options] -doSaf .
while read population bamfiles; do
    realSFS <Output_Path>/${population}.saf.idx > <Output_Path>/${population}.sfs;
done < tmp_${SLURM_ARRAY_TASK_ID} # Each array reads one line of the input file.

rm tmp_${SLURM_ARRAY_TASK_ID} # Remove temporary files created during the analysis.
