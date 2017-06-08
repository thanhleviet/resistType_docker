#!/bin/bash
#kmercount in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified May 23, 2014"
	echo ""
	echo "Usage:	kmercount in=<input> out=<read output> hist=<histogram output>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file.  'out' and 'hist' are both optional."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in2=null		Second input file for paired reads"
	echo "extra=null		Additional files to use for input (generating hash table) but not for output"
	echo "fastareadlen=2^31	Break up FASTA reads longer than this.  Can be useful when processing scaffolded genomes"
	echo "tablereads=-1		Use at most this many reads when building the hashtable (-1 means all)"
	echo "kmersample=1		Process every nth kmer, and skip the rest"
	echo "readsample=1		Process every nth read, and skip the rest"
	echo ""
	echo "Output parameters:"
	echo "hist=null		Specify a file to output the depth histogram"
	echo "histlen=10000		Max depth displayed on histogram"
	echo "reads=-1		Only process this number of reads, then quit (-1 means all)"
	echo "sampleoutput=true	Use sampling on output as well as input (not used if sample rates are 1)"
	echo "printcoverage=false	Only print coverage information instead of reads"
	echo "useheader=false	Append coverage info to the read's header"
	echo "minmedian=0		Don't output reads with median coverage below this"
	echo "minaverage=0		Don't output reads with average coverage below this"
	echo "zerobin=false		Set to true if you want kmers with a count of 0 to go in the 0 bin instead of the 1 bin in histograms."
	echo "				Default is false, to prevent confusion about how there can be 0-count kmers."
	echo "				The reason is that based on the 'minq' and 'minprob' settings, some kmers may be excluded from the bloom filter."
	echo ""
	echo "Hashing parameters:"
	echo "k=31			Kmer length (values under 32 are most efficient, but arbitrarily high values are supported)"
	echo "cbits=8			Bits per cell in bloom filter; must be 2, 4, 8, 16, or 32.  Maximum kmer depth recorded is 2^cbits."
	echo " 			Large values decrease accuracy for a fixed amount of memory."
	echo "hashes=4		Number of times a kmer is hashed.  Higher is slower."
	echo "  			Higher is MORE accurate if there is enough memory, and LESS accurate if there is not enough memory."
	echo "prefilter=false	True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable."
	echo "prehashes=2		Number of hashes for prefilter."
	echo "passes=1		More passes can sometimes increase accuracy by iteratively removing low-depth kmers"
	echo "minq=7			Ignore kmers containing bases with quality below this"
	echo "minprob=0.5		Ignore kmers with overall probability of correctness below this"
	echo "threads=X		Spawn exactly X hashing threads (default is number of logical processors).  Total active threads may exceed X by up to 4."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

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

kmercount() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.KmerCoverage prefilter=true bits=16 interleaved=false $@"
	echo $CMD >&2
	$CMD
}

kmercount "$@"
