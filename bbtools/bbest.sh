#!/bin/bash
#bbest in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Calculates EST (expressed sequence tags) capture by an assembly from a sam file.
Designed to use BBMap output generated with these flags: k=13 maxindel=100000 customtag ordered

Usage: bbest.sh in=<sam file> out=<stats file>

Parameters:
in=<file>             Specify a sam file (or stdin) containing mapped ests.
out=<file>            Specify the output stats file (default is stdout).
ref=<file>            Specify the reference file (optional).
est=<file>            Specify the est fasta file (optional).
fraction=<0.98>       Min fraction of bases mapped to ref to be 
                      considered 'all mapped'.

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

function bbest() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load samtools
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.SamToEst $@"
#	echo $CMD >&2
	$CMD
}

bbest "$@"
