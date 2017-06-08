#!/bin/bash
#filterbyname in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified December 5, 2014

Description:  Filters reads by name.

Usage:  filterbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.
Other parameters and their defaults:

include=f            Set to 'true' to include the filtered names rather than excluding them.
ow=t                (overwrite) Overwrites files that already exist.
app=f               (append) Append to files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
names=              A list of strings or files.  The files can have one name per line, or
                    be a standard read file (fasta, fastq, or sam).

Java Parameters:
-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'
To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx800m"
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
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function filterbyname() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP driver.FilterReadsByName $@"
	echo $CMD >&2
	$CMD
}

filterbyname "$@"
