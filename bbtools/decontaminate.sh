#!/bin/bash
#decontaminate in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell."
	echo "Last modified January 2, 2015"
	echo ""
	echo "Description:  Decontaminates multiplexed assemblies via normalization and mapping."
	echo ""
	echo "Usage:  decontaminate.sh reads=<file,file> ref=<file,file> out=<directory>"
	echo ""
	echo "Input Parameters:"
	echo "reads=<file,file>     Input reads; required parameter."
	echo "ref=<file,file>       Input assemblies; required parameter."
	echo "interleaved=auto      True forces paired/interleaved input; false forces single-ended mapping."
	echo "                      If not specified, interleaved status will be autodetected from read names."
	echo "unpigz=t              Spawn a pigz (parallel gzip) process for faster decompression.  Requires pigz to be installed."
	echo "touppercase=t         (tuc) Convert lowercase letters in reads to upper case (otherwise they will not match the reference)."
	echo ""
	echo "Output Parameters:"
	echo "pigz=f                Spawn a pigz (parallel gzip) process for faster compression.  Requires pigz to be installed."
	echo "tmpdir=.              Write temp files here.  By default is uses the system's $TMPDIR or current directory."
	echo "outdir=.              Write ouput files here."
	echo ""
	echo "Mapping Parameters:"
	echo "kfilter=55            Set to a positive number N to require minimum N contiguous matches for a mapped read."
	echo ""
	echo "Filtering Parameters:"
	echo "minc=3.5              Min average coverage to retain scaffold."
	echo "minp=20               Min percent coverage to retain scaffold."
	echo "minr=18               Min mapped reads to retain scaffold."
	echo "minl=500              Min length to retain scaffold."
	echo "ratio=1.2             Contigs will not be removed by minc unless the coverage changed by at least this factor.  0 disables this filter."
	echo "mapraw=t              Set true to map the unnormalized reads.  Required to filter by 'ratio'."
	echo ""
	echo "Normalization Parameters:"
	echo "mindepth=2            Min depth of reads to keep."
	echo "target=20             Target normalization depth."
	echo "hashes=4              Number of hashes in Bloom filter."
	echo "passes=1              Normalization passes."
	echo "minprob=0.5           Min probability of correctness to add a kmer."
	echo "depthpercentile=0.75  (dp) Percentile to use for depth proxy (0.5 means median)."
	echo "ecc=f                 Error-correction."
	echo "aecc=f                Agressive error-correction."
	echo "cecc=f                Conservative error-correction."
	echo "prefilter=f           Prefilter, for large datasets."
	echo "filterbits=32         (fbits) Bits per cell in primary filter."
	echo "prefilterbits=2       (pbits) Bits per cell in prefilter."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx                  This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "                      -Xmx20g will specify 20 gigs of RAM, and -Xmx800m will specify 800 megs.  The max is typically 85% of physical memory."
	echo ""
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"


decontaminate() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z $z2 -cp $CP jgi.DecontaminateByNormalization $@"
	echo $CMD >&2
	$CMD
}

decontaminate "$@"
