#!/bin/bash
#estherfilter.sh <query> <reference> <cutoff>"

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  BLASTs queries against reference, and filters out hits with
              scores less than 'cutoff'.  The score is taken from column 12
              of the BLAST output.  The specific BLAST command is:
              blastall -p blastn -i QUERY -d REFERENCE -e 0.00001 -m 8

Usage:  estherfilter.sh <query> <reference> <cutoff>

For example:

estherfilter.sh reads.fasta genes.fasta 1000 > results.txt

'fasta' can be used as a fourth argument to get output in Fasta format.  Requires more memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx3200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

estherfilter() {
	#module load oracle-jdk/1.7_64bit
	#module load blast
	local CMD="java $EA $z -cp $CP driver.EstherFilter $@"
	echo $CMD >&2
	$CMD
}

estherfilter "$@"
