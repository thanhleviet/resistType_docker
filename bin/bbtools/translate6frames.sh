#!/bin/bash
#translate6frames in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified July 29, 2014"
	echo ""
	echo "Description:  Translates nucleotide sequences to all 6 amino acid frames."
	echo ""
	echo "Usage:	translate6frames.sh in=<input file> out=<output file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Main input. in=stdin.fq will pipe from stdin."
	echo "in2=<file>       	Input for 2nd read of pairs in a different file."
	echo "interleaved=auto 	(int) t/f overrides interleaved autodetection."
	echo "qin=auto         	Input quality offset: 33 (Sanger), 64, or auto."
	echo "reads=-1         	If positive, quit after processing X reads or pairs."
	echo ""
	echo "Output parameters:"
	echo "out=<file>       	(outnonmatch) Write reads here that do not contain kmers matching the database.  'out=stdout.fq' will pipe to standard out."
	echo "out2=<file>      	(outnonmatch2) Use this to write 2nd read of pairs to a different file."
	echo "overwrite=t      	(ow) Grant permission to overwrite files."
	echo "ziplevel=2       	(zl) Compression level; 1 (min) through 9 (max)."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	Output quality offset: 33 (Sanger), 64, or auto."
	echo "tag=t            	Tag read id with the frame, adding e.g. ' fr1'"
	echo "frames=6         	Only print this many frames.  e.g. if you already know the sense, set 'frames=3'"
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

z="-Xmx2g"
z2="-Xms2g"
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
	freeRam 2000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

translate6frames() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP jgi.TranslateSixFrames $@"
	echo $CMD >&2
	$CMD
}

translate6frames "$@"
