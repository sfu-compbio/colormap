#!/bin/bash

# Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)

if [ "$#" -lt 5 ]; then
   echo "USAGE: ./runCorr.sh pacbio.fasta illumina.fastq outDirectory outPrefix threads"
   exit
fi

EXEDIR=$(dirname "${BASH_SOURCE[0]}")
CURRDIR=`pwd`
LRFILE=$(realpath "${1}")
SRFILE=$(realpath "${2}")
OUTDIR=$(realpath "${3}")
OUTPRE=$(basename "${4}")
NTHREAD=$5
BWAOPTIONS="-aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 16 -w 40 -r 1 -D 0 -y 20 -L 30,30 -T 2.5 -t $NTHREAD"

BWA=$EXEDIR/bin/bwa-proovread
SAM=$EXEDIR/bin/samtools
SPALG=$EXEDIR/bin/spCorrection
OEAALG=$EXEDIR/bin/oeaCorrection

UNCORR="$OUTDIR/${OUTPRE}_uncorr"
ITER1="$OUTDIR/${OUTPRE}_iter1"
ITER2="$OUTDIR/${OUTPRE}_iter2"

mkdir -p "$OUTDIR"
# rm -f "$3"/*
ln -s -f "$LRFILE" "${UNCORR}.fasta"

"$BWA" index "${UNCORR}.fasta"
"$BWA" mem $BWAOPTIONS "${UNCORR}.fasta" "$SRFILE" > "${ITER1}.sam"
rm "${UNCORR}.fasta.amb" "${UNCORR}.fasta.ann" "${UNCORR}.fasta.bwt" "${UNCORR}.fasta.pac" "${UNCORR}.fasta.sa"
"$SAM" view -bS "${ITER1}.sam" | "$SAM" sort -@ $NTHREAD - | "$SAM" view -h -o "${ITER1}.sort.sam" -
rm "${ITER1}.sam"
"$SPALG" -l "${UNCORR}.fasta" -a "${ITER1}.sort.sam" -t $NTHREAD > "${ITER1}.fasta"
rm "${ITER1}.sort.sam"

"$BWA" index "${ITER1}.fasta"
"$BWA" mem $BWAOPTIONS "${ITER1}.fasta" "$SRFILE" > "${ITER2}.sam"
rm "${ITER1}.fasta.amb" "${ITER1}.fasta.ann" "${ITER1}.fasta.bwt" "${ITER1}.fasta.pac" "${ITER1}.fasta.sa"
"$SAM" view -bS "${ITER2}.sam" | "$SAM" sort -@ $NTHREAD - | "$SAM" view -h -o "${ITER2}.sort.sam" -
rm "${ITER2}.sam"
"$SPALG" -l "${ITER1}.fasta" -a "${ITER2}.sort.sam" -t $NTHREAD > "${ITER2}.fasta"
rm "${ITER2}.sort.sam"
