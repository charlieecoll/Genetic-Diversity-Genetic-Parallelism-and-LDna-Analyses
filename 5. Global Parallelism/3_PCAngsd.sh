#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to calculate a PCA for each of the clusters of global parallelism with PCAngsd.
## The script is automated in arrays so each array will calculate the PCA for one of the clusters. It
## expects a file in which each row contains information on a cluster: the first column is the path
## and file to the beagle file of the cluster of parallelism and the second column the name of the cluster.

sed -n ${SLURM_ARRAY_TASK_ID}p <Path_and_File_to_be_Read_by_Arrays> > tmp_${SLURM_ARRAY_TASK_ID} # Set temporary files for the array.

while read file_path name; do
    python <Path>/pcangsd.py -beagle ${file_path} -o <Output_Path>/${name}_PCA -threads 64
done < tmp_${SLURM_ARRAY_TASK_ID} # Each array reads one line of the input file.

rm tmp_${SLURM_ARRAY_TASK_ID} # Remove temporary files created during the analysis.
