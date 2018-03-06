#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "USAGE: ./colormap_sge.sh CONFIG"
    exit
fi

job_id="1"

DIR_COLORMAP=`cat $1 | grep "CONF_COLORMAP" | awk '{print $2}'`
SP_ALG="${DIR_COLORMAP}/runCorr2.sh"
CHUNK_FILE=`cat $1 | grep "CONF_CHUNK" | awk '{print $2}'`
SR_FILE=`cat $1 | grep "CONF_SRFILE" | awk '{print $2}'`
NUM_THREAD=`cat $1 | grep "CONF_THREAD" | awk '{print $2}'`

LR_FILE=`cat ${CHUNK_FILE} | head -n ${job_id} | tail -n 1`

"${SP_ALG}" "${LR_FILE}" "${SR_FILE}" ${NUM_THREAD} 1
