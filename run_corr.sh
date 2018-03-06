#!/bin/bash

# Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)

if [ "$#" -lt 3 ]; then
    echo "USAGE: ./runCorr.sh long.fasta short.fastq threads [numIter(1)]"
    exit
fi

getFullDir_file(){
    echo "$(cd "$(dirname "$1")"; pwd)"
}
getFullDir_dir(){
    echo "$(mkdir -p "$1"; cd "$1"; pwd)"
}
getFullFile(){
    echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}
getFullBase(){
    tmpFull=`getFullFile "$1"`
    echo "${tmpFull%.*}"
}
getBase(){
    tmpFull=`getFullFile "$1"`
    tmpName=$(basename "$tmpFull")
    echo "${tmpName%.*}"
}

DIR_SCRIPT=`getFullDir_file "${BASH_SOURCE[0]}"`
FILE_LR=`getFullFile "$1"`
FILE_SR=`getFullFile "$2"`
BASE_LR=`getFullBase "$1"`
DIR_OUT="${BASE_LR}_tmp"
NAME_LR=`getBase "$1"`
NUM_THREAD="$3"
NUM_ITER="1"
if [ "$#" -gt 3 ]; then
    NUM_ITER="$4"
fi
if [ "${NUM_ITER}" -lt 1 ]; then
    NUM_ITER="1"
fi

TMP_UNCORR="${FILE_LR}"
BWAOPTIONS="-aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 16 -w 40 -r 1 -D 0 -y 20 -L 30,30 -T 2.5 -t ${NUM_THREAD}"
BWA="${DIR_SCRIPT}/bin/bwa-proovread"
SAMTOOLS="${DIR_SCRIPT}/bin/samtools"
SPALG="${DIR_SCRIPT}/bin/spCorrection"

# make the directory
mkdir -p "${DIR_OUT}"

# perform ${NUM_ITER} iterations of correction
for ii in `seq 1 "${NUM_ITER}"`; do
    ln -s -f "${TMP_UNCORR}" "${DIR_OUT}/${NAME_LR}_ITER${ii}.uncorr.fasta"
    TMP_UNCORR="${DIR_OUT}/${NAME_LR}_ITER${ii}.uncorr.fasta"
    # mapping
    "${BWA}" index "${TMP_UNCORR}"
    TMP_SAM="${DIR_OUT}/${NAME_LR}_ITER${ii}.sam"
    "${BWA}" mem ${BWAOPTIONS} "${TMP_UNCORR}" "${FILE_SR}" > "${TMP_SAM}"
    rm "${TMP_UNCORR}.amb" "${TMP_UNCORR}.ann" "${TMP_UNCORR}.bwt" "${TMP_UNCORR}.pac" "${TMP_UNCORR}.sa"
    echo "[NOTE] sorting SAM file: ${TMP_SAM}"
    # sorting
    "${SAMTOOLS}" view -bS -@ ${NUM_THREAD} "${TMP_SAM}" | "${SAMTOOLS}" sort -@ ${NUM_THREAD} - | "${SAMTOOLS}" view -h -@ ${NUM_THREAD} -o "${TMP_SAM}.sorted.sam" -
    rm "${TMP_SAM}"
    # correcting
    if [ "${ii}" -lt "${NUM_ITER}" ]; then
        TMP_CORR="${DIR_OUT}/${NAME_LR}_ITER${ii}.corr.fasta"
        echo "[NOTE] correcting long reads: ${TMP_UNCORR}"
        "${SPALG}" -l "${TMP_UNCORR}" -a "${TMP_SAM}.sorted.sam" -t ${NUM_THREAD} > "${TMP_CORR}"
        TMP_UNCORR="${TMP_CORR}"
    fi
done
# final correction
echo "[NOTE] correcting long reads: ${FILE_LR}"
"${SPALG}" -l "${TMP_UNCORR}" -a "${TMP_SAM}.sorted.sam" -t ${NUM_THREAD} > "${BASE_LR}.corr.fasta"
rm -rf ${DIR_OUT}
