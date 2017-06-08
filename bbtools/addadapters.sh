#!/bin/bash
#reformat in=<infile> out=<outfile>

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 1, 2014"
	echo ""
	echo "Description:  Randomly adds adapters to a file, or grades a trimmed file."
	echo ""
	echo "Usage:  addadapters.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> adapters=<file>"
	echo ""
	echo "in2 and out2 are for paired reads and are optional."
	echo "If input is paired and there is only one output file, it will be written interleaved."
	echo "Other parameters and their defaults:"
	echo ""
	echo "ow=f         		(overwrite) Overwrites files that already exist."
	echo "int=f		 	(interleaved) Determines whether INPUT file is considered interleaved."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo "add         		Add adapters to input files.  Default mode."
	echo "grade         		Evaluate trimmed input files."
	echo "adapters=<file>   	Fasta file of adapter sequences."
	echo "literal=<sequence>	Comma-delimited list of adapter sequences."
	echo "left         		Adapters are on the left (3') end of the read."
	echo "right         		Adapters are on the right (5') end of the read.  Default mode."
	echo "adderrors=t  		Add errors to adapters based on the quality scores."
	echo "addpaired=t  		Add adapters to the same location for read 1 and read 2."
	echo "arc=f        		Add reverse-complemented adapters as well as forward."
	echo "rate=0.5     		Add adapters to this fraction of reads."
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

function rename() {
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.AddAdapters $@"
	echo $CMD >&2
	$CMD
}

rename "$@"
