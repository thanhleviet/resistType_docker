#!/bin/bash
#matrixtocolumns in1=<infile> in2=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Transforms two matched identity matrices into 2-column format,
              one row per entry, one column per matrix.

Usage:  matrixtocolumns.sh in1=<matrix1> in2=<matrix2> out=<file>

Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

matrixtocolumns() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP driver.CorrelateIdentity $@"
	echo $CMD >&2
	$CMD
}

matrixtocolumns "$@"
