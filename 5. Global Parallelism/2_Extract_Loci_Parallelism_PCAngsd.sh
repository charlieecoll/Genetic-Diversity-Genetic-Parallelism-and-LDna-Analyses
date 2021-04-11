#!/bin/bash -l

## Carla Coll Costa, C3.

## Script to extract the loci of the clusters of global parallelism from the beagle file.
## The script creates one new beagle file for each cluster of global parallelism: each cluster (or new
## beagle file) contains a subset of loci of the beagle file with all the variants. The script is automated
## in arrays so each array will create a beagle file for one of the clusters. It expects a file with different
## rows in which each row has information on one cluster: the first column contains the path and name of a .txt file
## with the loci information on that cluster and the second column the name of the cluster. The file with the loci
## information consists on a .txt file where each line has the chromosome and position of a SNP in the following
## format: chromosome_position.

sed -n ${SLURM_ARRAY_TASK_ID}p <Path_and_File_to_be_Read_by_Arrays> > tmp_${SLURM_ARRAY_TASK_ID} # Set temporary files for the array.

# Generate a beagle file for each cluster
while read file_path cluster; do
    while read position; do
        grep -w -m1 $position <Beagle_File_with_All_Variants> >> <Output_Path_of_Cluster_Files>/${cluster}.beagle; # Create a beagle file for each of the clusters
    done < ${file_path}
done < tmp_${SLURM_ARRAY_TASK_ID} # Each array reads one line of the input file.

rm tmp_${SLURM_ARRAY_TASK_ID} # Remove temporary files created during the analysis.

# Get the heather of the general beagle file
head -n 1 <Beagle_File_with_All_Variants> > <Output_File_with_Heather>

# Merge "header" and content of each beagle file for each of the clusters, and compress them to final .gz file
for i in `seq 1 12`; do
    cat <Output_File_with_Heather> <Output_Path_of_Cluster_Files>/Cluster_${i}.beagle | gzip > <Output_Path_of_Final_Cluster_Files>/Cluster_${i}.beagle.gz;
done
