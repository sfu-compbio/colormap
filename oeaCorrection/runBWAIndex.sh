#!/bin/bash

# Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)

EXEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -lt 1 ]; then
   echo "USAGE: runBWAIndex.sh ref.fa"
   exit
fi

BWA=$EXEDIR/bwa-proovread

$BWA index $1
