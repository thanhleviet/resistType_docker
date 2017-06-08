#!/bin/bash
#mapnt in=<infile> out=<outfile>

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 16, 2014"
	echo "This script requires at least -Xmx106GB (120G node on Genepool)."
	echo ""
	echo "Description:  Maps sequences to the nt database."
	echo ""
	echo "Usage:	mapnt.sh in=<input file> outu=<clean output file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a fasta, fastq, or sam file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "interleaved=auto 	(int) If true, forces fastq input to be paired and interleaved."
	echo "trim=t           	Trim read ends to remove bases with quality below minq."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "untrim=t           	Undo the trimming after mapping."
	echo "trimq=10           	Trim quality threshold."
	echo "ziplevel=4       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "outm=<file>       	File to output the reads that mapped to human."
	echo ""
	echo "***** All BBMap parameters can be used; run bbmap.sh for more details. *****"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

function mapnt() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA -Xmx106g -Xms106g -Xmn4g -cp $CP align2.BBMap ow maxindel=80 k=11 minratio=0.45 usemodulo path=/global/projectb/sandbox/gaag/bbtools/nt/ build=2 pigz=f unpigz=f zl=4 qtrim=rl trimq=10 untrim idtag printunmappedcount local noheadersequences ambig=all $@"
	echo $CMD >&2
	$CMD
}

mapnt "$@"
