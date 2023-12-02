#!/bin/bash
set -ueox pipefail
# Log Location on Server.
# Assign date to variable.
today=$(date +"%Y-%m-%d_%H%M")


FIRSTDIR=''
SECONDDIR=''
THIRDDIR=''
THREADS=6

while getopts hr:t:o:a: flag
do
    case "${flag}" in
        r) OUT=${OPTARG};;
        o) FIRSTDIR=${OPTARG};;
        t) SECONDDIR=${OPTARG};;
        a) THIRDDIR=${OPTARG};;
        h)
                echo "Used for the help menu"
                echo "Usage: deeptools_bw.sh -l <dir to export log file> -r <dir of bam file> -o <outputdir> -p num.threads"
                echo "Made to be used with GNU parallel at this point"
                echo "Default log location is [ $LOG_LOCATION ]"
                exit
                ;;
    esac
done


exec 2>&1

echo "$OUT"
echo "$FIRSTDIR"
echo "$SECONDDIR"
echo "$THIRDDIR"

