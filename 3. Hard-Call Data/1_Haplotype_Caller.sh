#!/bin/bash -l

# Paolo Momigliano.

## Script to hard-call the variants from the .bam data.
## The script is automated so each array will hard-call the data for each .bam file. It expects a file
## in which each line correspond to a sample. In the first column there is the name of the sample and in
## the second column the path to the .bam file of that sample.

sed -n ${SLURM_ARRAY_TASK_ID}p <Path_and_File_to_be_Read_by_Arrays> > <Output_Path>/tmp_${SLURM_ARRAY_TASK_ID} # Set temporary files for the array.

time=`date +%F-%T`

module load biokit # Load module needed.
module load gatk # Load module needed.
echo -e \#\!\/\bin\/\bash > <Output_Path>/10_Join_CombineGVCFs.sh
echo reference=<Reference_Genome_with_Path> >> <Output_Path>/10_Join_CombineGVCFs.sh
echo -e gatk4 CombineGVCFs \\ >> <Output_Path>/10_Join_CombineGVCFs.sh
echo -e -R <Reference_Genome_with_Path> \\ >> <Output_Path>/10_Join_CombineGVCFs.sh


while read sample_name path; do \
   gatk --java-options "-Xmx32G"  HaplotypeCaller \
        -R <Reference_Genome_with_Path> \
        -I $path \
        -O <Output_Path>/${sample_name}.gvcf \
        -ERC GVCF \
        --native-pair-hmm-threads 16 \
        -ploidy 2
        
        echo -e $sample_name"\t"$gvcf >> <Output_Path>/gvcf.sample_map
        echo -e --variant $sample_name.gvcf \\ >> <Output_Path>/10_Join_CombineGVCFs.sh
        
done < <Output_Path>/tmp_${SLURM_ARRAY_TASK_ID} # Each array reads one line of the input file.

rm -f <Output_Path>/tmp_${SLURM_ARRAY_TASK_ID} # Remove temporary files created during the analysis.