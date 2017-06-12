#!/bin/bash

f=$1
if [ ! -s ${f}_1.fastq.gz ]
then
    fastq-dump --origfmt -I --split-files $f    --gzip 
fi

