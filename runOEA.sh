#!/bin/bash

# Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)

if [ "$#" -lt 5 ]; then
   echo "USAGE: ./runOEA.sh pacbio.fasta illumina.fastq outDirectory outPrefix threads"
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
# BWAOPTIONS="-aY -t $NTHREAD"

BWA=$EXEDIR/bin/bwa-proovread
SAM=$EXEDIR/bin/samtools
SPALG=$EXEDIR/bin/spCorrection
OEAALG=$EXEDIR/bin/oeaCorrection

CORR="$OUTDIR/${OUTPRE}_corr"
OEA="$OUTDIR/${OUTPRE}_oea"

mkdir -p "$OUTDIR"
# # rm -f $3/*
ln -s -f "$LRFILE" "${CORR}.fasta"

$BWA index "${CORR}.fasta"
$BWA mem $BWAOPTIONS "${CORR}.fasta" "$SRFILE" > "${OEA}.sam"
rm "${CORR}.fasta.amb" "${CORR}.fasta.ann" "${CORR}.fasta.bwt" "${CORR}.fasta.pac" "${CORR}.fasta.sa"
$SAM view -bS "${OEA}.sam" | $SAM sort -@ $NTHREAD - | $SAM view -h -o "${OEA}.sort.sam" -
rm "${OEA}.sam"
$OEAALG -s "$SRFILE" -l "${CORR}.fasta" -a "${OEA}.sort.sam" -t $NTHREAD > "${OEA}.fasta"
rm "${OEA}.sort.sam"
