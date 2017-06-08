#!/bin/bash
#readlength in=<infile>

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015
Description:  Generates a length histogram of input reads.

Usage:	readlength.sh in=<input file>

in=<file>    	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in.
in2=<file>   	Use this if 2nd read of pairs are in a different file.
out=<file>   	Write the histogram to this file.  Default is stdout.
bin=10       	Set the histogram bin size.
max=80000    	Set the max read length to track.
round=f      	Places reads in the closest bin, rather than the highest bin of at least readlength.
nzo=f        	(nonzeroonly) Do not print empty bins.
reads=-1     	If nonnegative, stop after this many reads.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx120m"
z2="-Xms120m"
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

stats() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.MakeLengthHistogram $@"
#	echo $CMD >&2
	$CMD
}

stats "$@"
