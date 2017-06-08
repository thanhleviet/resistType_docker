#!/bin/bash
#bbfakereads in=<infile> out=<outfile>

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 5, 2014"
	echo ""
	echo "Description:  Generates fake read pairs from ends of contigs or single reads."
	echo ""
	echo "Usage:  bbfakereads.sh in=<file> out=<outfile> out2=<outfile2>"
	echo ""
	echo "Out2 is optional; if there is only one output file, it will be written interleaved."
	echo "Other parameters and their defaults:"
	echo ""
	echo "ow=f         		(overwrite) Overwrites files that already exist."
	echo "zl=4            	(ziplevel) Set compression level, 1 (low) to 9 (max)."
	echo "fastawrap=100    	Length of lines in fasta output."
	echo "tuc=f    		(touppercase) Change lowercase letters in reads to uppercase."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo "qfin=<.qual file>	Read qualities from this qual file, for the reads coming from 'in=<fasta file>'"
	echo "qfout=<.qual file>	Write qualities from this qual file, for the reads going to 'out=<fasta file>'"
	echo "qfout2=<.qual file>	Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'"
	echo "verifyinterleaved=f	(vint) When true, checks a file to see if the names look paired.  Prints an error message if not."
	echo "tossbrokenreads=f	(tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash."
	echo ""
	echo "Faking parameters:"
	echo "length=250     	Generate reads of this length."
	echo "minlength=1       	Don't generate reads shorter than this."
	echo "overlap=0       	If you set overlap, then reads will by variable length, overlapping by 'overlap' in the middle."
	echo "identifier=null       (id) Output read names are prefixed with this."
	echo "addspace=t            Set to false to omit the  space before /1 and /2 of paired reads."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Supported input formats are fastq, fasta, fast+qual, scarf, and bread (BBMap's native format)"
	echo "Supported output formats are fastq, fasta, fast+qual, bread"
	echo "Supported compression formats are gz, zip, and bz2"
	echo "To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'"
	echo "To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx600m"
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

function bbfakereads() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP jgi.FakeReads $@"
	echo $CMD >&2
	$CMD
}

bbfakereads "$@"
