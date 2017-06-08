#!/bin/bash
#dedupe in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell and Jonathan Rood"
	echo "Last modified December 12, 2014"
	echo ""
	echo "Description:  Accepts one or more files containing sets of sequences (reads or scaffolds)."
	echo "Removes duplicate sequences, which may be specified to be exact matches, subsequences, or sequences within some percent identity."
	echo "Can also find overlapping sequences and group them into clusters."
	echo ""
	echo "Usage:	dedupe.sh in=<file or stdin> out=<file or stdout>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file.  With no output parameter, data will be written to stdout."
	echo "If 'out=null', there will be no output, but statistics will still be printed."
	echo "You can also use 'dedupe <infile> <outfile>' without the 'in=' and 'out='."
	echo ""
	echo "I/O Parameters"
	echo ""
	echo "in=<file,file>   	A single file or a comma-delimited list of files."
	echo "out=<file>       	Destination for all output contigs."
	echo "pattern=<file>   	Clusters will be written to individual files, where the '%' symbol in the pattern is replaced by cluster number."
	echo "outd=<file>      	Optional; removed duplicates will go here."
	echo "csf=<file>       	(clusterstatsfile) Write a list of cluster names and sizes."
	echo "dot=<file>       	(graph) Write a graph in dot format.  Requires 'fo' and 'pc' flags."
	echo "threads=auto          (t) Set number of threads to use; default is number of logical processors."
	echo "overwrite=t           (ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "showspeed=t           (ss) Set to 'f' to suppress display of processing speed."
	echo "minscaf=0             (ms) Ignore contigs/scaffolds shorter than this."
	echo "interleaved=auto      If true, forces fastq input to be paired and interleaved."
	echo "ziplevel=2            Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo ""
	echo "Output Format Parameters"
	echo ""
	echo "storename=t           (sn) Store scaffold names (set false to save memory)."
#	echo "addpairnum=f          Add .1 and .2 to numeric id of read1 and read2."
	echo "storequality=t        (sq) Store quality values for fastq assemblies (set false to save memory)."
	echo "uniquenames=t         (un) Ensure all output scaffolds have unique names.  Uses more memory."
	echo "numbergraphnodes=t    (ngn) Label dot graph nodes with read numbers rather than read names."
	echo "sort=f                Sort output by scaffold length (otherwise it will be random)."
	echo "                      'a' for ascending, 'd' for descending, 'f' for false (no sorting)."
	echo "renameclusters=f      (rnc) Rename contigs to indicate which cluster they are in."
	echo "printlengthinedges=f  (ple) Print the length of contigs in edges."
	echo ""
	echo "Processing Parameters"
	echo ""
	echo "absorbrc=t            (arc) Absorb reverse-complements as well as normal orientation."
	echo "absorbmatch=t         (am) Absorb exact matches of contigs."
	echo "absorbcontainment=t   (ac) Absorb full containments of contigs."
#	echo "absorboverlap=f       (ao) Absorb (merge) non-contained overlaps of contigs (TODO)."
	echo "findoverlap=f         (fo) Find overlaps between contigs (containments and non-containments).  Necessary for clustering."
	echo "uniqueonly=f          (uo) If true, all copies of duplicate reads will be discarded, rather than keeping 1."
	echo "rmn=f                 (requirematchingnames) If true, both names and sequence must match."
	echo "usejni=f              (jni) Do alignments in C code, which is faster, if an edit distance is allowed."
	echo "                      This will require compiling the C code; details are in /jni/README.txt."
	echo ""
	echo "Subset Parameters"
	echo ""
	echo "subsetcount=1         (sstc) Number of subsets used to process the data; higher uses less memory."
	echo "subset=0              (sst) Only process reads whose ((ID%subsetcount)==subset)."
	echo ""
	echo "Clustering Parameters"
	echo ""
	echo "cluster=f             (c) Group overlapping contigs into clusters."
	echo "pto=f                 (preventtransitiveoverlaps) Do not look for new edges between nodes in the same cluster."
	echo "minclustersize=1      (mcs) Do not output clusters smaller than this."
	echo ""
	echo "Cluster Postprocessing Parameters"
	echo ""
	echo "processclusters=f     (pc) Run the cluster processing phase, which performs the selected operations in this category."
	echo "fixmultijoins=t       (fmj) Remove redundant overlaps between the same two contigs."
	echo "removecycles=t        (rc) Remove all cycles so clusters form trees."
	echo "cc=t                  (canonicizeclusters) Flip contigs so clusters have a single orientation."
	echo "fcc=f                 (fixcanoncontradictions) Truncate graph at nodes with canonization disputes."
	echo "foc=f                 (fixoffsetcontradictions) Truncate graph at nodes with offset disputes."
	echo "mst=f                 (maxspanningtree) Remove cyclic edges, leaving only the longest edges that form a tree."
	echo ""
	echo "Overlap Detection Parameters"
	echo ""
	echo "exact=t               (ex) Only allow exact symbol matches.  When false, an 'N' will match any symbol."
	echo "touppercase=t         (tuc) Convert input bases to upper-case; otherwise, lower-case will not match."
	echo "maxsubs=0             (s) Allow up to this many mismatches (substitutions only, no indels).  May be set higher than maxedits."
	echo "maxedits=0            (e) Allow up to this many edits (subs or indels).  Higher is slower."
	echo "minidentity=100       (mid) Absorb contained sequences with percent identity of at least this (includes indels)."
	echo "minlengthpercent=0    (mlp) Smaller contig must be at least this percent of larger contig's length to be absorbed."
	echo "minoverlappercent=0   (mop) Overlap must be at least this percent of smaller contig's length to cluster and merge."
	echo "minoverlap=200        (mo) Overlap must be at least this long to cluster and merge."
	echo "depthratio=0          (dr) When non-zero, overlaps will only be formed between reads with a depth ratio of at most this."
	echo "                      Should be above 1.  Depth is determined by parsing the read names; this information can be added"
	echo "                      by running KmerNormalize (khist.sh, bbnorm.sh, or ecc.sh) with the flag 'rename'"
	echo "k=31                  Seed length used for finding containments and overlaps.  Anything shorter than k will not be found."
	echo "numaffixmaps=1        (nam) Number of prefixes/suffixes to index per contig. Higher is more sensitive, if edits are allowed."
#	echo "ignoreaffix1=f        (ia1) Ignore first affix (for testing)."
#	echo "storesuffix=f         (ss) Store suffix as well as prefix.  Automatically set to true when doing inexact matches."
	echo ""
	echo "Other Parameters"
	echo ""
	echo "forcetrimleft=-1      (ftl) If positive, trim bases to the left of this position (exclusive, 0-based)."
	echo "forcetrimright=-1     (ftr) If positive, trim bases to the right of this position (exclusive, 0-based)."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo "-Xms       		If you use the -Xmx flag, also set -Xms to the same value."
	echo ""
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

dedupe() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z $z2 -cp $CP jgi.Dedupe $@"
	echo $CMD >&2
	$CMD
}

dedupe "$@"
