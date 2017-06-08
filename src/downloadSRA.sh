#!/bin/bash
#$ -cwd
#$ -V

f=$1
if [ ! -s ${f}_1.fq ]
then
    fastq-dump --origfmt -I --split-files $f    --gzip 
fi

(mkdir -p fastq/$f
mv ${f}_1.fastq.gz fastq/$f/reads1.fq.gz
mv ${f}_2.fastq.gz fastq/$f/reads2.fq.gz) &
