#!/bin/bash
#seal in=<file> out=<file> ref=<ref file>

usage(){
echo "
Written by Brian Bushnell
Last modified January 13, 2015

Description:  Performs high-speed alignment-free sequence quantification,
by counting the number of long kmers that match between a read and
a set of reference sequences.  Designed for RNA-seq with alternative splicing.

Usage:  seal.sh in=<input file> ref=<file,file,file...> rpkm=<file>

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

out=<file>          (outmatch) Write reads here that contain kmers matching
                    the reference. 'out=stdout.fq' will pipe to standard out.
out2=<file>         (outmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outu=<file>         (outunmatched) Write reads here that do not contain kmers 
                    matching the database.
outu2=<file>        (outunmatched2) Use this to write 2nd read of pairs to a 
                    different file.
stats=<file>        Write statistics about which contamininants were detected.
refstats=<file>     Write statistics on a per-reference-file basis.
rpkm=<file>         Write RPKM for each reference sequence (for RNA-seq).
dump=<file>         Dump kmer tables to a file, in fasta format.
nzo=t               Only write statistics about ref sequences with nonzero hits.
overwrite=t         (ow) Grant permission to overwrite files.
showspeed=t         (ss) 'f' suppresses display of processing speed.
ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
fastawrap=80        Length of lines in fasta output.
qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
statscolumns=5      (cols) Number of columns for stats output, 3 or 5.
                    5 includes base counts.
rename=f            Rename reads to indicate which sequences they matched.
refnames=f          Use names of reference files rather than scaffold IDs.
trd=f               Truncate read and ref names at the first whitespace.
ordered=f           Set to true to output reads in same order as input.

Processing parameters:

k=31                Kmer length used for finding contaminants.  Contaminants 
                    shorter than k will not be found.  k must be at least 1.
rcomp=t             Look for reverse-complements of kmers in addition to 
                    forward kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to 
                    increase sensitivity in the presence of errors.
minkmerhits=1       (mkh) Reads with more than this many contaminant kmers 
                    will be discarded.
hammingdistance=0   (hdist) Maximum Hamming distance from ref kmers (subs only).
                    Memory use is proportional to (3*K)^hdist.
editdistance=0      (edist) Maximum edit distance from ref kmers (subs and 
                    indels).  Memory use is proportional to (8*K)^edist.
forbidn=f           (fn) Forbids matching of read kmers containing N.  
                    By default, these will match a reference 'A' if hdist>0
                    or edist>0, to increase sensitivity.
match=all           Determines when to quit looking for kmer matches.  Values:
                         all:    Attempt to match all kmers in each read.
                         first:  Quit aftet the first matching kmer.
                         unique: Quit after the first uniquely matching kmer.
ambiguous=random    (ambig) Set behavior on ambiguously-mapped reads (with an
                    equal number of kmer matches to multiple sequences).
                         first:  Use the first best-matching sequence.
                         toss:   Consider unmapped.
                         random: Select one best-matching sequence randomly.
                         all:    Use all best-matching sequences.

Speed and Memory parameters:

threads=auto        (t) Set number of threads to use; default is number of 
                    logical processors.
prealloc=f          Preallocate memory in table.  Allows faster table loading 
                    and more efficient memory usage, for a large reference.
rskip=1             Skip reference kmers to reduce memory usage.
                    1 means use all, 2 means use every other kmer, etc.
qskip=1             Skip query kmers to increase speed.  1 means use all.
speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
                    reads and reference.  Increases speed and reduces memory.
Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.

Trimming/Masking parameters:

qtrim=f             Trim read ends to remove bases with quality below minq.
                    Performed AFTER looking for kmers.  Values: 
                         t (trim both ends), 
                         f (neither end), 
                         r (right end only), 
                         l (left end only).
trimq=6             Regions with average quality BELOW this will be trimmed.
minlength=1         (ml) Reads shorter than this after trimming will be 
                    discarded.  Pairs will be discarded only if both are shorter.
maxlength=          Reads longer than this after trimming will be discarded.
                    Pairs will be discarded only if both are longer.
minavgquality=0     (maq) Reads with average quality (after trimming) below 
                    this will be discarded.
maxns=-1            If non-negative, reads with more Ns than this 
                    (after trimming) will be discarded.
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

Java Parameters:

-Xmx                This will be passed to Java to set memory usage, overriding 
                    the program's automatic memory detection. -Xmx20g will specify 
                    20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"	
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

seal() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z $z2 -cp $CP jgi.Seal $@"
	echo $CMD >&2
	$CMD
}

seal "$@"
