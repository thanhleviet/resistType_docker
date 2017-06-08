#!/bin/bash
#bbwrap in=<infile> out=<outfile>

usage(){
	echo "BBWrap v34.x"
	echo "Last modified September 18, 2014"
	echo ""
	echo "Description:  Wrapper for BBMap to allow multiple input and output files for the same reference."
	echo ""
	echo "To index:  	bbwrap.sh ref=<reference fasta>"
	echo "To map:    	bbwrap.sh in=<file,file,...> out=<file,file,...>"
	echo "To map without an index:  	bbwrap.sh ref=<reference fasta> in=<file,file,...> out=<file,file,...> nodisk"
	echo ""
	echo "This will probably not work with stdin and stdout, or histogram output."
	echo ""
	echo "Other Parameters:"
	echo ""
	echo "in=<file>        	Input sequences to mask. 'in=stdin.fa' will pipe from standard in."
	echo "mapper=bbmap     	Select mapper.  May be BBMap, BBMapPacBio, or BBMapPacBioSkimmer."
	echo "append=f         	Append to files rather than overwriting them.  If append is enabled,"
	echo "                 	multiple input files can write to a single output file, if there is exactly one output file."
	echo ""
	echo "***** All BBMap parameters can be used; run bbmap.sh for more details. *****"
	echo ""
	echo "Java Parameters:"
	echo "-Xmx             		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "                 		-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "This list is not complete.  For more information, please consult $DIR""docs/readme.txt"
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
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

bbwrap() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP align2.BBWrap build=1 overwrite=true fastareadlen=500 $@"
	echo $CMD >&2
	$CMD
}

bbwrap "$@"
