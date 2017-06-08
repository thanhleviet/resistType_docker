#!/bin/bash
#demuxbyname in=<infile> out=<outfile>

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified November 3, 2014"
	echo ""
	echo "Description:  Demultiplexes reads based on their name (suffix or prefix) into multiple files."
	echo ""
	echo "Usage:  demuxbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string...>"
	echo ""
	echo "in2 and out2 are for paired reads and are optional."
	echo "If input is paired and there is only one output file, it will be written interleaved."
	echo "Output filenames MUST contain a '%' symbol."
	echo "Other parameters and their defaults:"
	echo ""
	echo "names=        		List of strings (or files containing strings) to parse from read names."
	echo "prefix=t      		Match prefix of read name.  If false, match suffix of read name."
	echo "ow=f         		(overwrite) Overwrites files that already exist."
	echo "app=f        		(append) Append to files that already exist."
	echo "zl=4            	(ziplevel) Set compression level, 1 (low) to 9 (max)."
	echo "int=f		 	(interleaved) Determines whether INPUT file is considered interleaved."
	echo "fastawrap=80    	Length of lines in fasta output."
	echo "fastareadlen=0   	Set to a non-zero number to break fasta files into reads of at most this length."
	echo "minscaf=1       	Ignore fasta reads shorter than this."
	echo "tuc=f    		(touppercase) Change lowercase letters in reads to uppercase."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo "qfake=30        	Quality value used for fasta to fastq demuxbynameting."
	echo "qfin=<.qual file>	Read qualities from this qual file, for the reads coming from 'in=<fasta file>'"
	echo "qfin2=<.qual file>	Read qualities from this qual file, for the reads coming from 'in2=<fasta file>'"
	echo "qfout=<.qual file>	Write qualities from this qual file, for the reads going to 'out=<fasta file>'"
	echo "qfout2=<.qual file>	Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'"
	echo "tossbrokenreads=f	(tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash."
	echo "ignorebadquality=f	(ibq) Fix out-of-range quality values instead of crashing with a warning."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Supported input formats are fastq, fasta, fast+qual, scarf, and bread (BBMap's native format)"
	echo "Supported output formats are fastq, fasta, fast+qual, bread, sam, and bam (bam only if samtools is installed)"
	echo "Supported compression formats are gz, zip, and bz2"
	echo "To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'"
	echo "To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx400m"
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

function demuxbyname() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP jgi.DemuxByName $@"
	echo $CMD >&2
	$CMD
}

demuxbyname "$@"
