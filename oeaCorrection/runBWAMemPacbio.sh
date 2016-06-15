#!/bin/bash

# Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)

EXEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -lt 3 ]; then
   echo "USAGE: runBWAMemPacbio.sh reads.fa ref.fa outName"
   exit
fi

BWA=$EXEDIR/bwa-proovread

# $BWA index $2
$BWA mem -x pacbio $2 $1 | awk '{if($3!="*") print;}' 1> $3
