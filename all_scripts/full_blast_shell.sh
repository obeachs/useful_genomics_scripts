#!/bin/bash
set -ueo pipefail
BLAST_SCRIPT='/Volumes/sesame/joerecovery/scripts/full_blast_analysis_script_full.py'
PARSE_SCRIPT='/Volumes/sesame/joerecovery/scripts/parsing_blast_analysis.py'
FASTA_SCRIPT='/Volumes/sesame/joerecovery/scripts/generate_fasta_for_alignments.py'

FASTA=$1
DB=${2:-/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_cdna_20101214_updated}
DIR=${3:-${PWD}}

echo 'Running BLAST analysis on: '
echo  "$FASTA"
echo 'against the database: '
echo "$DB"
echo 'And saving results in the directory: '
echo "$DIR"

#Setting up directories and filenames 
BASE_FA="${FASTA##*/}"
BASE_FA=$(echo "$BASE_FA" | cut -f 1 -d '.')
BASE_DB="${DB##*/}"
mkdir -p "$DIR"/"$BASE_FA"_BLAST_TO_"$BASE_DB"
cp $FASTA "$DIR"/"$BASE_FA"_BLAST_TO_"$BASE_DB"
cd $DIR

FILE_LOC=${FASTA%/*}

echo "BEGINNING INITIAL BLAST"
python3.10 "${BLAST_SCRIPT}" --cdnafile $FASTA --outfolder ${DIR}/"$BASE_FA"_BLAST_TO_"$BASE_DB" --blastdatabase ${DB}


cd "$DIR"/"$BASE_FA"_BLAST_TO_"$BASE_DB"

mkdir -p alignments
echo "Parsing the BLAST analysis"
python3.10 $PARSE_SCRIPT --queryfasta $FASTA --refcdna $DB --blastoutput *dedup*
cd  alignments
pwd
echo 'Generating fastas and alignment files'

python3.10 $FASTA_SCRIPT --fasta "$DIR"/"$BASE_FA"_BLAST_TO_"$BASE_DB"/*parsed*

for i in $(ls *matches.fa | cut -f1 -d'.'); 
    do clustalo --force  -i $i.fa -o $i'_clustalo.fa';
done

