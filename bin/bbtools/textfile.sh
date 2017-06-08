#!/bin/bash

function usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Displays contents of a text file.

Usage:  textfile.sh <file> <start line> <stop line>

Start line and stop line are zero-based.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

function tf() {
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA -Xmx120m -cp $CP fileIO.TextFile $@"
	echo $CMD >&2
	$CMD
}

tf "$@"
