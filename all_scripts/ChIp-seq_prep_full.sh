#!/usr/bin/env bash
#This is just so that the program stops if it hits an error
set -ueo pipefail
#Just want to very simply download the SRR reads from NCBI and then carry out
#fqc quality check
#REF=$1
DIRECTORY=$1

for trimmed in $(ls $DIRECTORY | grep trimmed)
do
  REF=~/refs/Hisat2/tair10
  R4=${trimmed%.*}
  hisat2 -p 10 --no-spliced-alignment -x ~/refs/Hisat2/tair10 -U $R4.fq -S $R4'_hisat2.sam'
done

for sam in $(ls $DIRECTORY | grep sam)
do
  SAM=${sam%.*}
  samtools view -Sb $SAM.sam > $SAM.bam
done

for bam in $(ls $DIRECTORY | grep .bam)
do
  BAM=${bam%.*}
  samtools sort $BAM.bam  $BAM"_sorted.bam"
done

for sam in $(ls $DIRECTORY | grep sam)
do
  SAM=${sam%.*}
  rm $SAM.sam
done

for sorted in $(ls $DIRECTORY | grep sorted.bam)
do
  R7=${sorted%.*}
  echo "Flagstatting"
  samtools flagstat $R7.bam $R7'.flagstat''
done


for sorted in $(ls $DIRECTORY | grep sorted.bam)
do
  R7=${sorted%.*}
  echo "# Running samtools rmdup to remove PCR duplicates from $sorted"
  samtools rmdup -s $R7.bam $R7'_rmdup.bam'
done


for rmdup in $(ls $DIRECTORY | grep rmdup)
do
  rmdup=${rmdup%.*}
  echo $rmdup beep bop
  bedtools genomecov -ibam $rmdup.bam -g ~/refs/TAIR10_Chr.all.fai -bg > $rmdup'_bedtools.bedgraph'
done

for sorted in $(ls $DIRECTORY | grep rmdup)
do
  BAM=${sorted%.*}
  echo $BAM
  samtools index $BAM.bam
done

#Needed some extra info here in order to make sure we won't get bam.bai.bam.bw etc
#for bai in $(ls $DIRECTORY | grep bai | cut -f1 -d'.')
#do
#  BAI=${bai%.*}
#  bamCoverage -b $BAI.bam -o $BAI'_bamCov.bw'
#done
