#!/usr/bin/env bash
#This is just so that the program stops if it hits an error
set -ueo pipefail
#Just want to very simply download the SRR reads from NCBI and then carry out
#fqc quality check
#REF=$1
READ=$1
READ2=$2



echo "READ1 = $READ"
ehco "READ2 = $READ2"

HI_REF=~/refs/Hisat2/tair10
R1=${READ%.*}
R2=${READ2%.*}

echo "Aligning $READ and $READ2 with hisat2"
  hisat2 -p 10 --no-spliced-alignment -x ~/refs/Hisat2/tair10 -1 $READ -2 $READ2 $R4'_hisat2.sam'
done

for trimmed in $(ls $DIRECTORY | grep trimmed | cut -f1 -d '_')
do
  REF=~/refs/TAIR10_Chr.all.fa
  R4=${trimmed}
  echo "Processing $R4"
  bwa mem -t 2 $REF $R4_1-trimmed.fq $R4_2-trimmed.fq> $R4'_bwa.sam'
done

#for trimmed in $(ls $DIRECTORY | grep .sam)
#  do
#   R4=${trimmed%.*}
#   echo "# Running samtools view on $trimmed to remove nonuniquely aligned reads"
#   samtools view -Sh $R4.sam|
#   grep -e "^@" -e "XM:i:[012][^0-9]" |
#   grep -v "XS:i:" > $R4'_filtered.sam'
#done

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
  echo "#Flagstatting"
  #rmdup -s for SE reads, -S for PE reads
  samtools flagstat $R7.bam $R7'.flagstat''
done


for sorted in $(ls $DIRECTORY | grep sorted.bam)
do
  R7=${sorted%.*}
  echo "# Running samtools rmdup to remove PCR duplicates from $sorted"
  #rmdup -s for SE reads, -S for PE reads
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
for bai in $(ls $DIRECTORY | grep bai | cut -f1 -d'.' )
do
  BAI=${bai%.*}
  bamCoverage -b $BAI.bam -o $BAI'_bamCov.bw'
done
