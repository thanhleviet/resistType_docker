#!/bin/bash
#repair in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Re-pairs reads that became disordered or had some mates eliminated.

Usage:  repair.sh in=<input file> out=<pair output> outs=<singleton output>

Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed.
Output may be stdout or a file.

Optional parameters (and their defaults)

in=<file>       The 'in=' flag is needed if the input file is not the first 
                parameter.  'in=stdin' will pipe from standard in.
in2=<file>      Use this if 2nd read of pairs are in a different file.
out=<file>      The 'out=' flag is needed if the output file is not the second
                parameter.  'out=stdout' will pipe to standard out.
out2=<file>     Use this to write 2nd read of pairs to a different file.
outs=<file>     (outsingle) Write singleton reads here.

overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.
fint=f          (fixinterleaving) Fixes corrupted interleaved files using read
                names.  Only use on files with broken interleaving - correctly
                interleaved files from which some reads were removed.
repair=t        (rp) Fixes arbitrarily corrupted paired reads by using read 
                names.  Uses much more memory than 'fint' mode.
ain=f           (allowidenticalnames) When detecting pair names, allows 
                identical names, instead of requiring /1 and /2 or 1: and 2:

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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

repair() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.SplitPairsAndSingles rp $@"
	echo $CMD >&2
	$CMD
}

repair "$@"
