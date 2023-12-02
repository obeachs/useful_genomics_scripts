#!/bin/bash
set -ueox pipefail
# Log Location on Server.
# Assign date to variable.
today=$(date +"%Y-%m-%d_%H%M")


FIRSTDIR=''
SECONDDIR=''
THIRDDIR=''
THREADS=6



while getopts ho:b:c:a: flag
do
    case "${flag}" in
        o) OUT=${OPTARG};;
        a) FIRSTDIR=${OPTARG};;
        h) 
        		echo "Used for the help menu"
                echo "Usage: multibigwigsummary_for_pca.sh -o FULL_FILE_PATH_AND_NAME -a DIRECTORIES_OF_BIGWIG_FILES"
                echo "This will produce an uncompressed file that can be used in R, but the # might have to be removed"
                exit
                ;;
        
    esac
done


exec 2>&1


#Needs this other output for some reason
multiBigWigSummary bins -b ${FIRSTDIR}/*bw --outRawCounts $OUT -o ${OUT}_binary


