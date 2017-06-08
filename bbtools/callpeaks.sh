#!/bin/bash
#callpeaks in=<infile>

usage(){
	echo "Calls peaks from a 2-column (x, y) histogram."
	echo "Written by Brian Bushnell"
	echo "Last modified October 28, 2014"
	echo ""
	echo "Usage:	callpeaks.sh in=<histogram file> out=<output file>"
	echo ""
	echo "in=<file>        	'in=stdin.fq' will pipe from standard in."
	echo "out=<file>       	Write the peaks to this file.  Default is stdout."
	echo "minHeight=2     	(h) Ignore peaks shorter than this."
	echo "minVolume=2     	(v) Ignore peaks with less area than this."
	echo "minWidth=2      	(w) Ignore peaks narrower than this."
	echo "minPeak=2       	(minp) Ignore peaks with an X-value below this."
	echo "maxPeak=BIG       	(maxp) Ignore peaks with an X-value above this."
	echo "maxPeakCount=8  	(maxpc) Print up to this many peaks (prioritizing height)."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
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

stats() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA -Xmx120m -cp $CP jgi.CallPeaks $@"
#	echo $CMD >&2
	$CMD
}

stats "$@"
