#!/bin/bash
#filterbarcodes in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified July 9, 2014"
	echo ""
	echo "Description: Filters barcodes by quality, and generates quality histograms."
	echo ""
	echo "Usage:	filterbarcodes.sh in=<file> out=<file> maq=<integer>"
	echo ""
	echo "Input may be stdin or a fasta or fastq file, raw or gzipped."
	echo "Output may be stdout or a file."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Reads that have already been muxed with barcode qualities using mergebarcodes.sh."
	echo "interleaved=auto 	(int) If true, forces fastq input to be paired and interleaved."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo ""
	echo "Output parameters:"
	echo "out=<file>       	Write filtered reads here.  'out=stdout.fq' will pipe to standard out."
	echo "cor=<file>       	Correlation between read and index qualities."
	echo "bqhist=<file>       	Barcode quality histogram by position."
	echo "baqhist=<file>      	Barcode average quality histogram."
	echo "bmqhist=<file>      	Barcode min quality histogram."
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo "maq=0            	Filter reads with barcode average quality less than this."
	echo "mmq=0            	Filter reads with barcode minimum quality less than this."
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

z="-Xmx200m"
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

filterbarcodes() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.CorrelateBarcodes $@"
	echo $CMD >&2
	$CMD
}

filterbarcodes "$@"
