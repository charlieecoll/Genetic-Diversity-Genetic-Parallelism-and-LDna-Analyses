#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to calculate a PCA from all the filtered hard-called data.

module load plink # Load modules needed.

export OMP_NUM_THREADS=1 # Change number of threads to do the calculations.

plink --vcf <Input_Path_and_File_Name_VCF_File> --allow-extra-chr --pca --out <Output_Path_and_File_Name_PCA> # Calculate eigenvectors and eigenvalues for PCA of all samples