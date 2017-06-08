#!/bin/bash
#bbduk in=<file> out=<file> ref=<ref file>

usage(){
echo "
Written by Brian Bushnell
Last modified January 13, 2015

Description:  Compares reads to the kmers in a reference dataset, optionally 
allowing an edit distance. Splits the reads into two outputs - those that 
match the reference, and those that don't. Can also trim (remove) the matching 
parts of the reads rather than binning the reads.

Usage:  bbduk.sh in=<input file> out=<output file> ref=<contaminant files>

Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed.
Output may be stdout or a file.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped 
fasta input, set in=stdin.fa.gz

Optional parameters (and their defaults)

Input parameters:

in=<file>           Main input. in=stdin.fq will pipe from stdin.
in2=<file>          Input for 2nd read of pairs in a different file.
ref=<file,file>     Comma-delimited list of reference files.
literal=<seq,seq>   Comma-delimited list of literal reference sequences.
touppercase=f       (tuc) Change all bases upper-case.
interleaved=auto    (int) t/f overrides interleaved autodetection.
qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
reads=-1            If positive, quit after processing X reads or pairs.

Output parameters:

out=<file>          (outnonmatch) Write reads here that do not contain 
                    kmers matching the database.  'out=stdout.fq' will pipe 
                    to standard out.
out2=<file>         (outnonmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outm=<file>         (outmatch) Write reads here that contain kmers matching
                    the database.
outm2=<file>        (outmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outs=<file>         (outsingle) Use this to write singleton reads whose mate 
                    was trimmed shorter than minlen.
stats=<file>        Write statistics about which contamininants were detected.
refstats=<file>     Write statistics on a per-reference-file basis.
rpkm=<file>         Write RPKM for each reference sequence (for RNA-seq).
dump=<file>         Dump kmer tables to a file, in fasta format.
duk=<file>          Write statistics in duk's format.
nzo=t               Only write statistics about ref sequences with nonzero hits.
overwrite=t         (ow) Grant permission to overwrite files.
showspeed=t         (ss) 'f' suppresses display of processing speed.
ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
fastawrap=80        Length of lines in fasta output.
qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
statscolumns=3      (cols) Number of columns for stats output, 3 or 5.
                    5 includes base counts.
rename=f            Rename reads to indicate which sequences they matched.
refnames=f          Use names of reference files rather than scaffold IDs.
trd=f               Truncate read and ref names at the first whitespace.
ordered=f           Set to true to output reads in same order as input.

Histogram output parameters:

bhist=<file>        Base composition histogram by position.
qhist=<file>        Quality histogram by position.
aqhist=<file>       Histogram of average read quality.
bqhist=<file>       Quality histogram designed for box plots.
lhist=<file>        Read length histogram.
gchist=<file>       Read GC content histogram.
gcbins=100          Number gchist bins.  Set to 'auto' to use read length.

Histograms for sam files only (requires sam format 1.4 or higher):

ehist=<file>        Errors-per-read histogram.
qahist=<file>       Quality accuracy histogram of error rates versus quality 
                    score.
indelhist=<file>    Indel length histogram.
mhist=<file>        Histogram of match, sub, del, and ins rates by read location.
idhist=<file>       Histogram of read count versus percent identity.
idbins=100          Number idhist bins.  Set to 'auto' to use read length.

Processing parameters:

k=27                Kmer length used for finding contaminants.  Contaminants 
                    shorter than k will not be found.  k must be at least 1.
rcomp=t             Look for reverse-complements of kmers in addition to 
                    forward kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to 
                    increase sensitivity in the presence of errors.
maxbadkmers=0       (mbk) Reads with more than this many contaminant kmers 
                    will be discarded.
hammingdistance=0   (hdist) Maximum Hamming distance from ref kmers (subs only).
                    Memory use is proportional to (3*K)^hdist.
editdistance=0      (edist) Maximum edit distance from ref kmers (subs 
                    and indels).  Memory use is proportional to (8*K)^edist.
hammingdistance2=0  (hdist2) Sets hdist for short kmers, when using mink.
editdistance2=0     (edist2) Sets edist for short kmers, when using mink.
forbidn=f           (fn) Forbids matching of read kmers containing N.
                    By default, these will match a reference 'A' if 
                    hdist>0 or edist>0, to increase sensitivity.
removeifeitherbad=t (rieb) Paired reads get sent to 'outmatch' if either is 
                    match (or either is trimmed shorter than minlen).  
                    Set to false to require both.
findbestmatch=f     (fbm) If multiple matches, associate read with sequence 
                    sharing most kmers.  Reduces speed.

Speed and Memory parameters:

threads=auto        (t) Set number of threads to use; default is number of 
                    logical processors.
prealloc=f          Preallocate memory in table.  Allows faster table loading 
                    and more efficient memory usage, for a large reference.
minrskip=1          (mns) Force minimal skip interval when indexing reference 
                    kmers.  1 means use all, 2 means use every other kmer, etc.
maxrskip=99         (mxs) Restrict maximal skip interval when indexing 
                    reference kmers. Normally all are used for scaffolds<100kb, 
                    but with longer scaffolds, up to K-1 are skipped.
rskip=              Set both minrskip and maxrskip to the same value.
                    If not set, rskip will vary based on sequence length.
qskip=1             Skip query kmers to increase speed.  1 means use all.
speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
                    reads and reference.  Increases speed and reduces memory.
Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.

Trimming/Filtering/Masking parameters:

ktrim=f             Trim reads to remove bases matching reference kmers.
                    Values: 
                            t (trim)
                            f (don't trim), 
                            r (trim right end), 
                            l (trim left end), 
                            n (convert to N instead of trimming).
                    Any non-whitespace character other than t, f, r, l, n: 
                    convert to that symbol rather than trimming, and process 
                    short kmers on both ends.
useshortkmers=f     (usk) Look for shorter kmers at read tips (only for 
                    k-trimming).  Enabling this will disable maskmiddle.
mink=6              Minimum length of short kmers.  Setting this automatically 
                    sets useshortkmers=t.
qtrim=f             Trim read ends to remove bases with quality below minq.
                    Performed AFTER looking for kmers.
                    Values: 
                            t (trim both ends), 
                            f (neither end), 
                            r (right end only), 
                            l (left end only).
trimq=6             Regions with average quality BELOW this will be trimmed.
minlength=10        (ml) Reads shorter than this after trimming will be 
                    discarded.  Pairs will be discarded if both are shorter.
maxlength=          Reads longer than this after trimming will be discarded.
                    Pairs will be discarded only if both are longer.
minavgquality=0     (maq) Reads with average quality (after trimming) below 
                    this will be discarded.
chastityfilter=f    (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
maxns=-1            If non-negative, reads with more Ns than this 
                    (after trimming) will be discarded.
otm=f               (outputtrimmedtomatch) Output reads trimmed to shorter 
                    than minlength to outm rather than discarding.
tp=0                (trimpad) Trim this much extra around matching kmers.
tbo=f               (trimbyoverlap) Trim adapters based on where paired 
                    reads overlap.
minoverlap=24       Require this many bases of overlap for overlap detection.
mininsert=25        Require insert size of at least this much for overlap 
                    detection (will automatically set minoverlap too).
tpe=f               (trimpairsevenly) When kmer right-trimming, trim both 
                    reads to the minimum length of either.
forcetrimleft=0     (ftl) If positive, trim bases to the left of this position
                    (exclusive, 0-based).
forcetrimright=0    (ftr) If positive, trim bases to the right of this position
                    (exclusive, 0-based).
forcetrimmod=0      (ftm) If positive, right-trim length to be equal to zero,
                    modulo this number.
restrictleft=0      If positive, only look for kmer matches in the 
                    leftmost X bases.
restrictright=0     If positive, only look for kmer matches in the 
                    rightmost X bases.

Entropy parameters:

entropy=-1          Set between 0 and 1 to filter reads with entropy below
                    that value.  Higher is more stringent.
entropywindow=50    Calculate entropy using a sliding window of this length.
entropyk=5          Calculate entropy using kmers of this length.

Java Parameters:

-Xmx                This will be passed to Java to set memory usage, overriding 
                    the program's automatic memory detection. -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.

There is a changelog at /bbmap/docs/changelog_bbduk.txt
Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"	
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
	freeRam 2000m 42
	z="-Xmx2G"
	z2="-Xms2G"
}
calcXmx "$@"

bbduk() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z $z2 -cp $CP jgi.BBDukF $@"
	echo $CMD >&2
	$CMD
}

bbduk "$@"
