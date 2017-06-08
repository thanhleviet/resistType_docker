#!/bin/bash

function usage(){
	echo "BBSplit / BBMap v34.x"
	echo "Written by Brian Bushnell, from Dec. 2010 - present"
	echo "Last modified December 11, 2014"
	echo ""
	echo "Description:  Maps reads to multiple references simultaneously."
	echo "Outputs reads to a file for the reference they best match, with multiple options for dealing with ambiguous mappings."
	echo ""
	echo "To index:  	bbsplit.sh build=<1> ref_x=<reference fasta> ref_y=<another reference fasta>"
	echo "To map:    	bbsplit.sh build=<1> in=<reads> out_x=<output file> out_y=<another output file>"
	echo ""
	echo "To be concise, and do everything in one command:"
	echo "bbsplit.sh ref=x.fa,y.fa in=reads.fq basename=o%.sam"
	echo "which is equivalent to"
	echo "bbsplit.sh build=1 in=reads.fq ref_x=x.fa ref_y=y.fa out_x=ox.sam out_y=oy.sam"
	echo ""
	echo "in=stdin will accept reads from standard in, and out=stdout will write to standard out,"
	echo "but file extensions are still needed to specify the format of the input and output files."
	echo "e.g. in=stdin.fa.gz will read gzipped fasta from standard in; out=stdout.sam.gz will write gzipped sam."
	echo ""	
	echo "Indexing Parameters (required when building the index):"
	echo "ref_<name>=<ref.fasta>  	Specify the reference sequence for the given name; e.g., ref_ecoli=ecoli.fasta"
	echo "                      	These can also be comma-delimited lists of files; e.g., ref_a=a1.fa,a2.fa,a3.fa"
	echo "build=<1>        		If multiple references are indexed in the same directory, each needs a unique build ID."
	echo ""
	echo "Input Parameters:"
	echo "build=<1>        		Designate index to use.  Corresponds to the number specified when building the index."
	echo "in=<reads.fq>    		Primary reads input; required parameter."
	echo "in2=<reads2.fq>  		For paired reads in two files."
	echo "qin=<auto>       		Set to 33 or 64 to specify input quality value ASCII offset."
	echo "interleaved=<auto>  		True forces paired/interleaved input; false forces single-ended mapping."
	echo "                 		If not specified, interleaved status will be autodetected from read names."
	echo ""
	echo "Mapping Parameters:"
	echo "maxindel=<20>    		Don't look for indels longer than this.  Lower is faster.  Set to >=100k for RNA-seq."
	echo "minratio=<0.9>   		Fraction of max alignment score required to keep a site.  Higher is faster."
	echo "minhits=<2>      		Minimum number of seed hits required for candidate sites.  Higher is faster."
	echo "ambiguous=<best> 		Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations)."
	echo "                 			best	(use the first best site)"
	echo "                 			toss	(consider unmapped)"
	echo "                 			random	(select one top-scoring site randomly)"
	echo "                 			all	(retain all top-scoring sites.  Does not work yet with SAM output)"
	echo "ambiguous2=<best> 		Set behavior only for reads that map ambiguously to multiple different references."
	echo "                 		Normal 'ambiguous=' controls behavior on all ambiguous reads;"
	echo "                 		Ambiguous2 excludes reads that map ambiguously within a single reference."
	echo "                 			best	(use the first best site)"
	echo "                 			toss	(consider unmapped)"
	echo "                 			all	(write a copy to the output for each reference to which it maps)"
	echo "                 			split	(write a copy to the AMBIGUOUS_ output for each reference to which it maps)"
	echo "trim=<true> 		 	Quality-trim ends to Q5 before mapping.  Options are 'l' (left), 'r' (right), and 'lr' (both)."
	echo "untrim=<true>  		Undo trimming after mapping.  Untrimmed bases will be soft-clipped in cigar strings."
	echo ""
	echo "Output Parameters:"
	echo "out_<name>=<file>		Output reads that map to the reference <name> to <file>."
	echo "basename=prefix%suffix	Equivalent to multiple out_%=prefix%suffix expressions, in which each % is replaced by the name of a reference file."
	echo "bs=<file>        		Write a shell script to 'file' that will turn the sam output into a sorted, indexed bam file."
	echo "scafstats=<file>     		Write statistics on how many reads mapped to which scaffold to this file."
	echo "refstats=<file>      		Write statistics on how many reads mapped to which reference to this file."
	echo ""
	echo "***** All BBMap parameters can be used; run bbmap.sh for more details. *****"
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
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

function bbsplit() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP align2.BBSplitter build=1 overwrite=true match=long fastareadlen=500 minhits=2 minratio=0.9 maxindel=20 trim=both untrim=true $@"
	echo $CMD >&2
	$CMD
}

bbsplit "$@"
