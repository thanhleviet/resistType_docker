#!/bin/bash
#mergeOTUs in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Merges coverage stats lines (from pileup) for the same OTU,
              according to some custom naming scheme.

Usage:  mergeOTUs.sh in=<file> out=<file>

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
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

function mergeOTUs() {
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA $z -cp $CP driver.MergeCoverageOTU $@"
	echo $CMD >&2
	$CMD
}

mergeOTUs "$@"
