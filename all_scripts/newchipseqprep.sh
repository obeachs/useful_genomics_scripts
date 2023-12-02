#!/bin/bash/

# initialize a variable with an intuitive name to store the name of the input fastq file


fq1=$1
fq2=$2
dir=${fq1%.*}
aligner=$3


genome=~/Documents/TAIR10_Chr_all



#Just extracting file names correctly and makinng directories for the intermediary files of the analysis
title1=`basename $fq1 .fq`
title2=`basename $fq2 .fq`
echo 'Analysing sample $title'
echo $dir/
echo $VARIABLE

mkdir -p $dir/fastqc_results/
mkdir -p $dir/bamfiles/


#Running fastQC on both of the files
echo 'Running fastQC on $title1 and $title2'
fastqc $fq1 $fq2 -o $dir/fastqc_results


# if [$aligner == 'bowtie']
# then
# 	#Using bowtie2 to align with high specificity
# 	bowtie2 -p 6 -q --local -x $genome -1 $fq1 -2 $fq2 -S $dir/bamfiles/$title1_$title2.sam
bowtie2 -p 6 -q --local -x $genome -1 $fq1 -2 $fq2 -S $dir/bamfiles/$title1_$title2.sam
#Using hisat2 to align with lower specificity, faster
#hisat2 -p 10 --no-spliced-alignment -x ~/refs/Hisat2/tair10 -1 $fq1 -2 $fq2 $dir/bamfiles/$title1_$title_hisat2.sam

#Samtools section
samtools view -h -S -b -@ 6 -o $dir/bamfiles/$title1_$title2_unsorted.bam $dir/bamfiles/$title1_$title2.sam
samtools sort $dir/bamfiles/$title1_$title2_unsorted.bam -o $dir/bamfiles/$title1_$title2_sorted.bam
samtools index $dir/bamfiles/$title1_$title2_sorted.bam

#Make bedgraph files for IGV visualisation
bedtools genomecov -ibam $dir/bamfiles/$title1_$title2_sorted.bam -g $genome -bg > $dir/bamfiles/$title1_$title2_sorted.begraph