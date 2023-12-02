#!/usr/bin/env bash
#This is just so that the program stops if it hits an error
set -ueo pipefail
#Just want to very simply download the SRR reads from NCBI and then carry out
#fqc quality check
#REF=$1
READ=$1
READ2=$2

REF=~/refs/Hisat2/tair10
R1=${READ%_*}

R2=${READ2%_*}

COMBINED = cut $READ -f1 -d'_'

echo $REF
echo $R1
echo $R2
echo $COMBINED
