#!/bin/bash
#merge in=<infile> out=<outfile>

function usage(){
echo "
BBMerge v5.5
Written by Brian Bushnell and Jonathan Rood
Last modified January 23, 2014

Description:  Merges paired reads into single reads by overlap detection.
With sufficient coverage, can also merge nonoverlapping reads using gapped kmers.

Usage for interleaved files:	bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
Usage for paired files:     	bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

Input may be stdin or a fasta, fastq, or scarf file, raw or gzipped.
Output may be stdout or a file.

Optional parameters (and their defaults)

Input parameters:
in=null              Primary input. 'in2' will specify a second file.
interleaved=auto     May be set to true or false to override autodetection of
                     whether the input file as interleaved.
reads=-1             Quit after this many read pairs (-1 means all).

Output parameters:
out=<file>           File for merged reads. 'out2' will specify a second file.
outu=<file>          File for unmerged reads. 'outu2' will specify a second file.
outinsert=<file>     File list of read names and their insert sizes.
hist=null            Insert length histogram output file.
nzo=t                Only print histogram bins with nonzero values.
ziplevel=2           Set to 1 (lowest) through 9 (max) to change compression
                     level; lower compression is faster.

Trimming/Filtering parameters:
qtrim=f              Trim read ends to remove bases with quality below minq.  
                     Performed BEFORE merging.
                     Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=10             Trim quality threshold.
minlength=0          (ml) Reads shorter than this after trimming, but before
                     merging, will be discarded. Pairs will be discarded only
                     if both are shorter.
maxlength=-1         Reads with longer insert sizes will be discarded.
trimonfailure=t      (tof) If detecting insert size by overlap fails,
                     the reads will be trimmed and this will be re-attempted.
tbo=f                (trimbyoverlap) Trim overlapping reads to remove 
                     rightmost (3') non-overlapping portion, instead of joining.
minavgquality=0      (maq) Reads with average quality below this (AFTER trimming)
                     will not be attempted to be merged.
maxexpectederrors=0  (mee) If positive, reads with more combined expected 
                     errors than this will not be attempted to be merged.

Processing Parameters:
join=t               Create merged reads.  If set to false, you can still 
                     generate an insert histogram.
useoverlap=t         Attempt merge based on paired read overlap.
minoverlap=12        Minimum number of overlapping bases to consider merging
                     reads.
minoverlap0=8        Minimum number of overlapping bases to consider for 
                     deciding ambiguity.
minoverlapinsert=25  Do not look for insert sizes shorter than this.
mininsert=35         Reads with insert sizes less than this (after merging) 
                     will be discarded.
minq=9               Ignore bases with quality below this.
maxq=41              Cap output quality scores at this.
margin=2             The best overlap must have at least 'margin' fewer 
                     mismatches than the second best.
mismatches=3         Allow up to this many mismatches in the overlap.
entropy=t            Increase the minimum overlap requirement for low-
                     complexity reads.
efilter=f            Ban overlaps with many more mismatches than expected.  
                     Default: true for loose/vloose, false otherwise.
kfilter=f            Ban overlaps that create novel kmers.  Does not seem to help.

Processing Modes (these are mutually exclusive macros that set other parameters):
strict=f             Decrease false positive rate and merging rate.
fast=f               Increase speed and slightly decrease merging rate.
normal=t             Default.
loose=f              Increase false positive rate and merging rate.
vloose=f             Greatly increase false positive rate and merging rate.
usejni=f             (jni) Do overlapping in C code, which is faster.  Requires
                     compiling the C code; details are in /jni/README.txt.

Java Parameters:
-Xmx                 This will be passed to Java to set memory usage, 
                     overriding the program's automatic memory detection.
                     For example, -Xmx400m will specify 400 MB RAM.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

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

function merge() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP jgi.BBMerge $@"
	echo $CMD >&2
	$CMD
}

merge "$@"
