#!/bin/bash

## Paolo Momigliano.

## Script to map the fq.gz files to the reference genome. 
## This script will also concatenate the data if we have more than one sequence file for individual
## (if individuals have been resequenced). 

## By pair of fastq (i.e. read groups) 
	## Step1: MAPPING EACH READ PAIR -> bwa mem
	## Step2a: FILL IN MATE COORDINATES AND INSERT SIZE FIELDS -> samtools fixmate -m
	## Step2b: SORT READS BY COORDINATES -> samtools sort (by coordinates)
	## Step2c: INDEX SAM/BAM FILES -> samtools index (tarvitaan seuraavaan)

## By sample
	## Step3: MERGE MULTIPLE FILES IN ONE OUTPUT -> samtools merge
	## Step4a: LIST OF REGIONS TO REALIGN -> gatk ReAlignerTargetCreator now set up to run with 30 cores, change -nt option if needed
	## Step4b: REALIGN AROUND INDELS -> gatk indelReAligner now set up to run with 30 cores, change -nt option if needed
	## Step4c: INDEX SAM/BAM FILES -> samtools index

debug=0

time=`date +%F-%T`
reference=/scratch/project_2003317/3s/1_Reference/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa #path to the reference genome. Remember to index prior to start!
path=/scratch/project_2003317/3s/4_Clone_Filter/ #path to the folder with fastq files

pop=$1 # population given as an argument on command line
sample=$2 # individual given as an argument on command line
BAMMAP=/scratch/project_2003317/3s/5_BAM/$pop/bamFileList3s ## tab-separated file to store paths of the output bamfiles

outdir=/scratch/project_2003317/3s/5_BAM/$pop/$sample # output directory
tmpdir=/scratch/project_2003317/3s/tmp/3s_map/$pop/$sample # temporal directory
sortdir=/scratch/project_2003317/3s/tmp/sort/ # dierctory for sorting

ls $path/$pop/$sample/*1.f*.gz 2>/dev/null


FILES1=( `ls $path/$pop/$sample/*1.f*.gz 2>/dev/null`" " )  # lists fastq files (mate) sequenced at BGI and MACROGEN
FILES2=( `ls $path/$pop/$sample/*2.f*.gz 2>/dev/null`" " ) # lists fastq files (pair) sequenced at BGI and MACROGEN
echo $FILES1
echo $FILES2

##	if false; then


mkdir -p $outdir
mkdir -p $tmpdir
	
cd $tmpdir
	
	
for ((i=0;i<${#FILES1[@]};++i)); do ## for loop over the read groups (each mate and pair is one readgroup)
    FLOWCELL=`zcat ${FILES1[i]} |head -1|cut -f3 -d":"`
    LANE=`zcat ${FILES1[i]} |head -1|cut -f4 -d":"`
    ID=$sample.$FLOWCELL.$LANE # genereate read group id
    
	echo $ID
    if true; then			
	
	logfile=$tmpdir/$sample""_analysis.log.$time
	
	out1=$tmpdir/$sample""_$ID""_aligned.bam # bwa mem out
	out2a=$tmpdir/$sample""_$ID""_fixmate.bam # fixmate out 
	out2b=$tmpdir/$sample""_$ID""_sorted.bam # sorted bam out
	out2c=$tmpdir/$sample""_$ID""_sorted.bai # indexed bam out
	out3=$tmpdir/$sample""_merged.bam # bam after merging
	out3b=$tmpdir/$sample""_merged.bai # index for merged bam
	out4a=$tmpdir/$sample""_target_intervals.list # targets for realignment
	out4b=$tmpdir/$sample"".bam # realigned bam
	out4c=$tmpdir/$sample"".bai # index for realigned bam

	
	echo $out1
	echo $out2a
	echo $out2b
	echo $out2c 
	# Step1: bwa
	
	date >> $logfile
	echo -e "\n3s: bwa\nRead group: "$ID >> $logfile
	echo -e "${FILES1[i]}\n${FILES2[i]}" >> $logfile
	
	if [ ! -e $out1 ] && [ ! -e $out2a ] && [ ! -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	
	  cat ${FILES1[i]} > file1.fq.gz
	  cat ${FILES2[i]} > file2.fq.gz
	  
	  ( bwa mem -t 1 -M -R "@RG\tID:$ID\tSM:$sample\tPL:ILLUMINA\tLB:LIB-1" \
		$reference \
		file1.fq.gz \
		file2.fq.gz  \
		| samtools view -h -b -o $out1 - ) 2>> $logfile
		
	fi
	
	
	
	# Step2a: samtools fixmate
	echo >> $logfile
	date >> $logfile
	echo -e "\n3s: samtools fixmate\nRead group: "$ID >> $logfile
	
	if [ -e $out1 ] && [ ! -e $out2a ] && [ ! -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	
	  samtools fixmate -m $out1 $out2a &>> $logfile  
	  if [ -e $out2a ] && [ $debug -eq 0 ]; then 
		  rm -f $out1 file1.fq.gz file2.fq.gz
	  fi
	fi
	
	
	# Step2b: samtools sort
	echo >> $logfile
	date >> $logfile
	echo -e "\n3s: samtools sort\nRead group: "$ID >> $logfile
		
	if [ -e $out2a ] && [ ! -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	  
	  samtools sort -@ 1 -T $sortdir/$sample""_$ID -O bam -o $out2b $out2a &>> $logfile
	  if [ -e $out2b ] && [ $debug -eq 0 ]; then 
		  rm -f $out2a 
	  fi
	fi
	
	# Step2c: samtools index
	echo >> $logfile
	date >> $logfile
	echo -e "\n3s: samtools index\nRead group: "$ID >> $logfile
	
	if [ -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	  samtools index  $out2b $out2c &>> $logfile
	fi 
    
    fi ## End of the debug statement
done ## end of the for loop over read groups

if true; then

debug=0

#By sample
#Step3: samtools merge
#Step4a: gatk ReAlignerTargetCreator
#Step4b:gatk indelReAligner
#Step5:samtools markdup
#Step5a: samtools index



# Step3 samtools merge
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools merge\n" >> $logfile
if [ ! -e $out3 ] && [ ! -e $out3b ] && [ ! -e $out4a ] && [ ! -e $out4b ]  ; then
  samtools merge $out3 $tmpdir/*_sorted.bam
  if [ -e $out3 ] && [ $debug -eq 0 ]; then
    rm -f $tmpdir/$sample""*_sorted.ba*
  fi
fi

# Step 3b samtools index
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools index\n" >> $logfile
if [ -e $out3 ] && [ ! -e $out3b ] && [ ! -e $out4a ] && [ ! -e $out4b ] ; then
    samtools index $out3 $out3b &>> $logfile
fi

# Step4a  gatk ReAlignerTargetCreator 
echo >> $logfile
date >> $logfile
echo -e "\n3s: gatk RealignerTargetCreator\n" >> $logfile

if [ -e $out3 ] && [ -e $out3b ] && [ ! -e $out4a ] && [ ! -e $out4b ]  ; then

  java -Xmx8G -jar /projappl/project_2003317/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R $reference \
	-I $out3 \
	-o $out4a \
	-nt 1 \
	&>> $logfile
	  
fi

# Step4: gatk indelReAligner
echo >> $logfile
date >> $logfile
echo -e "\n3s: gatk indelReAligner\n" >> $logfile

if [ -e $out3 ] && [ -e $out3b ] && [ -e $out4a ] && [ ! -e $out4b ] && [ ! -e $out4c ]; then

  java -Xmx8G -jar /projappl/project_2003317/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $reference \
	-I $out3 \
	-targetIntervals $out4a \
	-o $out4b \
	&>> $logfile


###	While debugging avoid removing all there
  # if [ -e $out3 ] && [ -e $out3b ] && [ -e $out4a ] && [ -e $out4b ] && [ $debug -eq 0 ]; then
	# rm -f $out3 $out3b $out4a
  # fi	
fi


#Step5:samtools markdup
# DUPLICATES ARE NOT REMOVED, THEY ARE ONLY MARKED!
# echo >> $logfile
# date >> $logfile
# echo -e "\n3s: samtools markdup\n" >> $logfile

# if [ -e $out4b ] && [ -e $out4c ] && [ ! -e $out4c] && [ ! -e $out5a ]; then

 #  samtools markdup $out4b $out5

# fi

#Step5a: samtools index
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools index\n" >> $logfile

if  [ -e $out4b ] && [ ! -e $out4c ]; then

    samtools index $out4b $out4c &>> $logfile
fi


##	below lines to copy the output and remove files. Leave out for now while troubleshootingrm
if [ -e $out4b ] && [ -e $out4c ] ; then
 echo -e "$sample\t$outdir/$sample"".bam" >> $BAMMAP
 mv $out4b $outdir
 mv $out4c $outdir
 rm -f $out3 $out3b $out4a
fi
echo >> $logfile
date >> $logfile
echo -e "\n3s: Pipeline finished\n" >> $logfile

mv $logfile $outdir
fi ## End of the debugging statement
	# fi

# Output: 1 bam and 1 bai file per sample!