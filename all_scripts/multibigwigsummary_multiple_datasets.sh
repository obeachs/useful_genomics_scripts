#!/bin/bash
set -ueox pipefail
# Log Location on Server.
# Assign date to variable.
today=$(date +"%Y-%m-%d_%H%M")


THREADS=6

while :; do
    case $1 in
        -a|--out) flag1            
        ;;
        -c|--dir1) optflag1="SET"            
        ;;
        -d|--dir2) optflag2="SET"            
        ;;
        -e|--dir3) optflag3="SET"            
        ;;
        *) break
    esac
    shift
done



echo "${flag1}"
echo "${optflag1}"
echo "${optflag2}"
echo "${optflag3}"


echo "DONE DONE DONE"