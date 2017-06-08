#!/bin/bash
#calctruequality in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified May 23, 2014"
	echo ""
	echo "Description:  Calculates the observed quality scores from a sam file."
	echo ""
	echo "Usage:	calctruequality.sh in=<file> out=<file> sam=<file,file,...file>"
	echo ""
	echo "Input may be stdin or a fasta or fastq file, raw or gzipped."
	echo "Output may be stdout or a file."
	echo "sam is optional, but may be a comma-delimited list of sam files to mask."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Sam file.  Must be in 1.4 format (with = and X cigar symbols, not M)."
	echo "reads=-1        	Stop after processing this many reads (if positive)."
	echo ""
	echo "Output parameters:"
	echo "overwrite=f      	(ow) Set to true to allow overwriting of existing files."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo ""
	echo "Other parameters:"
	echo "pigz=t    	  	Use pigz to compress.  If argument is a number, that will set the number of pigz threads."
	echo "unpigz=t    	  	Use pigz to decompress."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx400m"
z2="-Xms400m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
}
calcXmx "$@"

calctruequality() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.CalcTrueQuality $@"
	echo $CMD >&2
	$CMD
}

calctruequality "$@"
