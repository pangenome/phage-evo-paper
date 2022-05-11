#!/usr/bin/env bash

function fastq_count_bases
# This function counts the number of bases on a
# BGZIPPED fastq file 
{
    zcat $1 | awk 'BEGIN {sum=0;} {if(NR%4==2) {sum+=length($0);}} END {print sum;}'
}

abundance=$1    # abundance of the P1 passage
bp=$2           # Total bases of the P1 assembly
fastq_input=$3  # Input BGZIPPED fastq file 

echo "${abundance}" "$(fastq_count_bases $fastq_input)" "${bp}" \
    | awk '{print ($1*$2)/$3}' # prints an abundance relative to P1
