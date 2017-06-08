#!/bin/bash
#msa in=<file> out=<file>

usage(){
echo "
Written by Brian Bushnell
Last modified December 5, 2014

Description:  Aligns a query sequence to reference sequences.
Outputs the best matching position per reference sequence.
If there are multiple queries, only the best-matching query will be used.

Usage:	msa.sh in=<file> out=<file> query=<literal,literal,...>

Parameters:

in=<file>           File containing reads. in=stdin.fa will pipe from stdin.
out=<file>          Matrix output. out=stdout will pipe to stdout.
literal=            A sequence of bases to match, or a comma-delimited list.

Java Parameters:

-Xmx                This will be passed to Java to set memory usage, overriding 
                    the program's automatic memory detection. -Xmx20g will specify 
                    20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

msa() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.FindPrimers $@"
	echo $CMD >&2
	$CMD
}

msa "$@"
