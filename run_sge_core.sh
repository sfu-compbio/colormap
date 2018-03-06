#!/bin/bash

if [ "$#" -lt 5 ]; then
    echo "USAGE: ./colormap_core.sh long.fasta short.fastq chunkSize minLen numThread"
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
DIR_COLORMAP="${DIR_SCRIPT}/colormap"
FILE_SR=`getFullFile "$2"`
FILE_LR=`getFullFile "$1"`
NAME_LR=`getBase "$1"`
PREF_SPLIT="${NAME_LR}_chunk"
BASE_LR=`getFullBase "$1"`
DIR_OUT="${BASE_LR}_tmp"
SIZE_CHUNK="$3"
LEN_MIN="$4"
NUM_THREAD="$5"

# make output directory
mkdir -p ${DIR_OUT}
# split the long read file into smaller chunks
"${DIR_SCRIPT}/bin/fastUtils" split -f <( "${DIR_SCRIPT}/bin/fastUtils" format -f "${FILE_LR}" -m ${LEN_MIN} ) -s ${SIZE_CHUNK} -p "${DIR_OUT}/${PREF_SPLIT}"
ls -1 "${DIR_OUT}"/"${PREF_SPLIT}"* > "${DIR_OUT}"/allChunks.fofn
NUM_JOBS=`cat "${DIR_OUT}"/allChunks.fofn | wc -l`
# create CONFIG file
echo "CONF_COLORMAP ${DIR_SCRIPT}" > "${DIR_OUT}"/CONFIG
echo "CONF_CHUNK ${DIR_OUT}"/allChunks.fofn >> "${DIR_OUT}"/CONFIG
echo "CONF_SRFILE ${FILE_SR}" >> "${DIR_OUT}"/CONFIG
echo "CONF_THREAD ${NUM_THREAD}" >> "${DIR_OUT}"/CONFIG
