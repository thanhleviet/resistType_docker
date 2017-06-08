#!/bin/bash
#mapPacBio in=<infile> out=<outfile> ref=<reference>

usage(){
	bash "$DIR"bbmap.sh
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

mapPacBio() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP align2.BBMapPacBio build=1 overwrite=true minratio=0.40 fastareadlen=6000 ambiguous=best minscaf=100 startpad=10000 stoppad=10000 midpad=6000 $@"
	echo $CMD >&2
	$CMD
}

mapPacBio "$@"
