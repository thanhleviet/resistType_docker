#!/bin/bash
#makechimeras in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified October 20, 2014"
	echo ""
	echo "Description:  Makes chimeric PacBio reads from nonchimeric reads."
	echo ""
	echo "Usage:	makechimeras.sh in=<input> out=<output> readsout=<integer>"
	echo ""
	echo "Input may be stdin or a file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo ""
	echo "Input Parameters:"
	echo "in=<file>		The input file containing nonchimeric reads."
	echo "unpigz=<true>		Decompress with pigz for faster decompression."
	echo ""
	echo "Output Parameters:"
	echo "out=<file>		Fasta output destination."
	echo "chimeras=-1		Number of chimeras to create."
	echo "forcelength=0		If a positive number X, one parent will be length X, and the other will be length-X."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
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
	freeRam 800m 82
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

makechimeras() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.MakeChimeras $@"
	echo $CMD >&2
	$CMD
}

makechimeras "$@"
