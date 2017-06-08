#!/bin/bash
#reformat in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified January 23, 2015

Description:  Reformats reads to change ASCII quality encoding, interleaving, file format, or compression format.
Optionally performs additional functions such as quality trimming, subsetting, and subsampling.
Supports sam, fastq, fasta, fasta+qual, scarf, gzip, zip.

Usage:  reformat.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.
Other parameters and their defaults:

ow=f                    (overwrite) Overwrites files that already exist.
app=f                   (append) Append to files that already exist.
zl=4                    (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f                   (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=80            Length of lines in fasta output.
fastareadlen=0          Set to a non-zero number to break fasta files into reads of at most this length.
minscaf=1               Ignore fasta reads shorter than this.
tuc=f                   (touppercase) Change lowercase letters in reads to uppercase.
qin=auto                ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto               ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
qfake=30                Quality value used for fasta to fastq reformatting.
qfin=<.qual file>       Read qualities from this qual file, for the reads coming from 'in=<fasta file>'
qfin2=<.qual file>      Read qualities from this qual file, for the reads coming from 'in2=<fasta file>'
qfout=<.qual file>      Write qualities from this qual file, for the reads going to 'out=<fasta file>'
qfout2=<.qual file>     Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'
verifypaired=f          (vpair) When true, checks reads to see if the names look paired.  Prints an error message if not.
verifyinterleaved=f     (vint) sets 'vpair' to true and 'interleaved' to true.
allowidenticalnames=f   (ain) When verifying pair names, allows identical names, instead of requiring /1 and /2 or 1: and 2:
tossbrokenreads=f       (tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash.
ignorebadquality=f      (ibq) Fix out-of-range quality values instead of crashing with a warning.
addslash=f              Append ' /1' and ' /2' to read names, if not already present.  Also add 'int=t' if the reads are interleaved.
rcomp=f                 (rc) Reverse-compliment reads.
rcompmate=f             (rcm) Reverse-compliment read 2 only.
mingc=0                 Discard reads with GC content below this.
maxgc=1                 Discard reads with GC content above this.
mappedonly=f            Toss unmapped reads.
unmappedonly=f          Toss mapped reads.
changequality=t         (cq) N bases always get a quality of 0 and ACGT bases get a min quality of 2.
fixquality=f            Quality scores above 41 are capped at 41.

Histogram output parameters:
bhist=<file>            Base composition histogram by position.
qhist=<file>            Quality histogram by position.
aqhist=<file>           Histogram of average read quality.
bqhist=<file>           Quality histogram designed for box plots.
lhist=<file>            Read length histogram.
gchist=<file>           Read GC content histogram.
gcbins=100              Number gchist bins.  Set to 'auto' to use read length.

Histograms for sam files only (requires sam format 1.4 or higher):
ehist=<file>            Errors-per-read histogram.
qahist=<file>           Quality accuracy histogram of error rates versus quality score.
indelhist=<file>        Indel length histogram.
mhist=<file>            Histogram of match, sub, del, and ins rates by read location.
idhist=<file>           Histogram of read count versus percent identity.
idbins=100              Number idhist bins.  Set to 'auto' to use read length.

Sampling parameters:
reads=-1                Set to a positive number to only process this many INPUT reads (or pairs), then quit.
samplerate=1            Randomly output only this fraction of reads; 1 means sampling is disabled.
sampleseed=-1           Set to a positive number to use that prng seed for sampling (allowing deterministic sampling).
samplereadstarget=0     (srt) Exact number of OUTPUT reads (or pairs) desired.
samplebasestarget=0     (sbt) Exact number of OUTPUT bases desired.
                        Important: srt/sbt flags should not be used with stdin, samplerate, qtrim, minlength, or minavgquality.

Trimming and filtering parameters:
qtrim=f                 Trim read ends to remove bases with quality below minq.
                        Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=6                 Regions with average quality BELOW this will be trimmed.
minlength=0             (ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.
maxlength=0             If nonzero, reads longer than this after trimming will be discarded.
breaklength=0           If nonzero, reads longer than this will be broken into multiple reads of this length.  Does not work for paired reads.
requirebothbad=t        (rbb) Only discard pairs if both reads are shorter than minlen.
minavgquality=0         (maq) Reads with average quality (after trimming) below this will be discarded.
chastityfilter=f        (cf) Reads with names  containing ' 1:Y:' or ' 2:Y:' will be discarded.
maxns=-1                If 0 or greater, reads with more Ns than this (after trimming) will be discarded.
forcetrimleft=0         (ftl) If nonzero, trim left bases of the read to this position (exclusive, 0-based).
forcetrimright=0        (ftr) If nonzero, trim right bases of the read after this position (exclusive, 0-based).
forcetrimmod=5          (ftm) If positive, trim length to be equal to zero modulo this number.
outsingle=<file>        (outs) If a read is longer than minlength and its mate is shorter, the longer one goes here.

Sam file reformatting options.  Note that most of these will require an indexed reference.
build=<integer>         Assign a genome's build id.  You can index like this: bbmap.sh ref=<file> build=1
sam=1.4                 Set to 1.4 to write Sam version 1.4 cigar strings, with = and X, or 1.3 to use M.
md=f                    Set to true to write MD tags.
xs=f                    Set to 'ss', 'fs', or 'us' to write XS tags for RNAseq using secondstrand, firststrand,
                        or unstranded libraries.  Needed by Cufflinks.  JGI mainly uses 'firststrand'.
stoptag=t               Set to true to write a tag indicating read stop location, prefixed by YS:i:
idtag=t                 Set to true to write a tag indicating percent identity, prefixed by YI:f:

Java Parameters:
-Xmx                    This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                        -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Supported input formats are fastq, fasta, fast+qual, scarf, and bread (BBMap's native format)
Supported output formats are fastq, fasta, fast+qual, bread, sam, and bam (bam only if samtools is installed)
Supported compression formats are gz, zip, and bz2
To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'
To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

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

function reformat() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP jgi.ReformatReads $@"
	echo $CMD >&2
	$CMD
}

reformat "$@"
