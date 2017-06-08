#!/bin/bash
#bbcountunique in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified September 19, 2014"
	echo ""
	echo "Description:  Generates a kmer uniqueness histogram, binned by file position."
	echo "There are 3 columns for single reads, 6 columns for paired:"
	echo "count      	number of reads or pairs processed"
	echo "r1_first   	percent unique 1st kmer of read 1"
	echo "r1_rand    	percent unique random kmer of read 1"
	echo "r2_first   	percent unique 1st kmer of read 2"
	echo "r2_rand    	percent unique random kmer of read 2"
	echo "pair       	percent unique concatenated kmer from read 1 and 2"
	echo ""
	echo ""
	echo "Usage:	bbcountunique.sh in=<input> out=<output>"
	echo ""
	echo "Input may be a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in2=null		Second input file for paired reads"
	echo "interleaved=auto	Set true/false to override autodetection of the input file as paired interleaved."
	echo "samplerate=1		Set to below 1 to sample a fraction of input reads."
	echo "reads=-1		Only process this number of reads, then quit (-1 means all)"
	echo ""
	echo "Output parameters:"
	echo "out=<file>        	File for output stats"
	echo ""
	echo "Processing parameters:"
	echo "k=20			Kmer length (range 1-31)."
	echo "interval=25000		Print one line to the histogram per this many reads."
	echo "cumulative=f		Show cumulative numbers rather than per-interval numbers."
	echo "percent=t		Show percentages of unique reads."
	echo "count=f			Show raw counts of unique reads."
	echo "printlastbin=f		(plb) Print a line for the final undersized bin."
	echo ""
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

z="-Xmx1g"
z2="-Xms1g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bbcountunique() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.CalcUniqueness $@"
	echo $CMD >&2
	$CMD
}

bbcountunique "$@"
