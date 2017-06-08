#!/bin/bash
#bbqc in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 5, 2014"
	echo ""
	echo "Description:  Performs quality-trimming; artifact, human, and phiX removal; adapter-trimming; error-correction and normalization."
	echo "Designed for Illumina fragment libraries only."
	echo ""
	echo "Usage:	bbqc.sh in=<input file> out=<output file>"
	echo ""
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Input reads."
	echo "in2=<file>       	Use this if 2nd read of pairs are in a different file."
	echo "ref=<file,file>  	Comma-delimited list of additional reference files for filtering."
	echo "artifactdb=<file>	Override default Illumina artifacts file."
	echo "rnadb=<file>		Override default rna spikein file."
	echo "dnadb=<file>		Override default dna spikein file."
	echo "phixref=<file>   	Override default phiX reference file."
	echo ""
	echo "Output parameters:"
	echo "path=null       		Set to the directory to use for all output files."
	echo "out=null			Read output file name.  If this is left blank, the input filename will be used with '.filtered' inserted before the extension."
	echo "scafstats=scaffoldStats.txt	Scaffold stats file name (how many reads matched which reference scaffold) ."
	echo "kmerstats=kmerStats.txt    	Kmer stats file name (duk-like output)."
	echo "log=status.log   		Progress log file name."
	echo "filelist=file-list.txt   	Progress log file name."
	echo "stats=filterStats.txt   	Overall stats file name."
	echo "reproduceName=reproduce.sh	Name of shellscript to reproduce these results."
	echo ""
	echo "Processing parameters:"
	echo "rna=f            	Set to 't' for rna libraries, 'f' for dna libraries."
	echo "phix=t           	Set to 't' to remove reads containing phiX kmers."
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "filterk=27       	Kmer length for finding contaminants.  Contaminants shorter than k will not be found.."
	echo "trimk=23         	Kmer length for linker/adapter trimming."
	echo "mink=11               Minimum kmer length for short kmers when trimming."
	echo "mapk=13           	Kmer length for mapping to human."
	echo "normalizek=27    	Kmer length for normalization/error correction."
	echo "rcomp=t          	Look for reverse-complements of kmers in addition to forward kmers."
	echo "maskmiddle=t     	(mm) Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors."
	echo "maxbadkmers=0    	(mbk) Reads with more than this many contaminant kmers will be discarded."
	echo "filterhdist=1    	Hamming distance used for filtering."
	echo "trimhdist=1      	Hamming distance used for trimming."
	echo "trimhdist2=      	Hamming distance used for trimming with short kmers.  If unset, trimhdist will be used."
	echo "mapindex=          	Remove contaminants by mapping to the index at this location; default is an HG19 index."
	echo "mapref=          	Remove contaminants by mapping to this fasta file.  Overrides mapindex flag."
	echo ""
	echo "Quality trimming parameters:"
	echo "qtrim=rl          	Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=10           	Trim quality threshold."
	echo "minlength=40     	(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter."
	echo "minlengthfraction=0.6	(mlf) Reads shorter than this fraction of original length after trimming will be discarded."
	echo "minavgquality=8  	(maq) Reads with average quality (before trimming) below this will be discarded."
	echo "maxns=1          	Reads with more Ns than this will be discarded."
	echo "forcetrimmod=5   	(ftm) If positive, trim length to be equal to zero modulo this number."
	echo ""
	echo "Normalization/Error Correction Paramters:"
	echo "ecc=f	         	Correct errors."
	echo "aecc=f	     		Aggressive error correction (fixes more errors; less conservative)."
	echo "cecc=f	     		Conservative error correction."
	echo "normalize=f           Normalize."
	echo "target=50	     	Target kmer depth after normalization."
	echo "min=-1	     		Throw away reads below this depth (default is 1/8th of target)."
	echo "max=-1	     		Retain all reads between above target but below this depth (default is 11/10ths of target)."
	echo "passes=2     		Override number of passes for normalization.  By default it is 2 for normalization and 1 for error-correction only."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo ""
	echo "*****   All additional parameters supported by BBDukF may also be used, and will be passed directly to BBDukF   *****"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
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
	freeRam 9500m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bbqc() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	export TZ="America/Los_Angeles" 
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP jgi.BBQC usejni $@"
	echo $CMD >&2
	$CMD
}

bbqc "$@"
