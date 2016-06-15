#!/bin/bash

# Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)

EXEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -lt 5 ]; then
   echo "USAGE: runMinia.sh workingDir reads.fasta kmerSize minAbundance prefix"
   exit
fi

cd $1
MINIA=$EXEDIR/minia
$MINIA -in $2 -kmer-size $3 -abundance-min $4 -out $5 -nb-cores 1
