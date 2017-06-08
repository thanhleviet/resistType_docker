BBMap readme by Brian Bushnell
Last updated January 23, 2015.
Please contact me at bbushnell@lbl.gov if you have any questions or encounter any errors.
BBMap and all other BBTools are free to use for noncommercial purposes, and investigators are free to publish results derived from them, as long as the source code is not published without explicit permission.
The BBTools package was written by Brian Bushnell, with the exception of the (optional, but faster) C, JNI, and MPI components, which were written by Jonathan Rood.

Special thanks for help with shellscripts goes to:
Alex Copeland (JGI), Douglas Jacobsen (JGI/NERSC), and sdriscoll (SeqAnswers).

This is the official release of BBMap, version 34.x


Basic Syntax:

(Using shellscript, on Genepool, which autodetects RAM to set -Xmx parameter.  You can also include a flag like '-Xmx31g' in the shellscript arguments to set RAM usage.)
To index:
bbmap.sh ref=<reference.fa>
To map:
bbmap.sh in=<reads.fq> out=<mapped.sam>

(without shellscript)
To index:
java -ea -Xmx31g -cp <PATH> align2.BBMap ref=<reference.fa>
To map:
java -ea -Xmx31g -cp <PATH> align2.BBMap in=<reads.fq> out=<mapped.sam>

...where "<PATH>" should indicate the path to the directory containing all the source code directories; e.g. "/usr/bin/bbmap/current"

Please note, the reference is only needed for building the index the first time; subsequently, just specify the build number which corresponds to that reference.
So for example the first time you map to e.coli you might specify "ref=ecoli_reference.fa build=3"; after that, just specify "build=3".
The index files would then be stored in ./ref/genome/3/ and ./ref/index/3/
Also, the -Xmx parameter should specify approximately 85% of the physical memory of the target machine; so, 21G for a 24GB node.  The process needs approximately 8 bytes per reference base (plus a several hundred MB overhead).


Advanced Syntax:


Indexing Parameters (required when building the index):
path=<.>        	Base directory to store index files.  Default is the local directory.  The index will always be placed in a subdirectory "ref".
ref=<ref.fasta> 	Use this file to build the index.  Needs to be specified only once; subsequently, the build number should be used.
build=<1>		Write the index to this location (build=1 would be stored in /ref/genome/1/ and /ref/index/1/).  Can be any integer.  This parameter defaults to 1, but using additional numbers allows multiple references to be indexed in the same directory.
k=<13>          	Use length 13 kmers for indexing.  Suggested values are 9-15, with lower typically being slower and more accurate.  13 is usually optimal.  14 is better for RNA-SEQ and very large references >4GB; 12 is better for PacBio and cross-species mapping.
midpad=<300>		Put this many "N" in between scaffolds when making the index.  300 is fine for metagenomes with millions of contigs; for a finished genome like human with 25 scaffolds, this should be set to 100000+ to prevent cross-scaffold mapping.
startpad=<8000> 	Put this many "N" at the beginning of a "chrom" file when making index.  It's best if this is longer than your longest expected read.
stoppad=<8000>		Put this many "N" at the end of a "chrom" file when making index.  It's best if this is longer than your longest expected read.
minscaf=<1>		Do not include scaffolds shorter than this when generating index.  Useful for assemblies with millions of fairly worthless unscaffolded contigs under 100bp.  There's no reason to make this shorter than the kmer length.
usemodulo=<f>		Throw away ~80% of kmers based on their remainder modulo a number.  Reduces memory usage by around 50%, and reduces sensitivity slightly.  Must be specified when indexing and when mapping.


Input Parameters:
path=<.>		Base directory to read index files.
build=<1>		Use the index at this location (same as when indexing).
in=<reads.fq>		Use this as the input file for reads.  Also accepts fasta.  "in=sequential length=200" will break a genome into 200bp pieces and map them to itself.  "in=stdin" will accept piped input.  The format of piped input can be specified with e.g. "in=stdin.fq.gz" or "in=stdin.fa"; default is uncompressed fastq.
in2=<reads2.fq> 	Run mapping paired, with reads2 in the file "reads2.fq"
			NOTE:  As a shorthand, "in=reads#.fq" is equivalent to "in=reads1.fq in2=reads2.fq"
interleaved=<auto>	Or "int". Set to "true" to run mapping paired, forcing the reads to be considered interleaved from a single input file.  By default the reader will try to determine whether a file is interleaved based on the read names; so if you don't want this, set interleaved=false.
qin=<auto>       	Set to 33 or 64 to specify input quality value ASCII offset.
fastareadlen=<500>	If fasta is used for input, breaks the fasta file up into reads of about this length.  Useful if you want to map one reference against another, since BBMap currently has internal buffers limited to 500bp.  I can change this easily if desired.
fastaminread=<1>	Ignore fasta reads shorter than this.  Useful if, say, you set fastareadlen=500, and get a length 518 read; this will be broken into a 500bp read and an 18bp read.  But it's not usually worth mapping the 18bp read, which will often be ambiguous.
maxlen=<0>        	Break long fastq reads into pieces of this length.
minlen=<0>       	Throw away remainder of read that is shorter than this.
fakequality=<-1>	Set to a positive number 1-50 to generate fake quality strings for fasta input reads.  Less than one turns this function off.
blacklist=<a.fa,b.fa>	Set a list of comma-delimited fasta files.  Any read mapped to a scaffold name in these files will be considered "blacklisted" and can be handled differently by using the "outm", "outb", and "outputblacklisted" flags.  The blacklist fasta files should also be merged with other fasta files to make a single combined fasta file; this combined file should be specified with the "ref=" flag when indexing.
touppercase=<f>		Set true to convert lowercase read bases to upper case.  This is required if any reads have lowercase letters (which real reads should never have).


Sampling Parameters:
reads=<-1>		Process at most N reads, then stop.  Useful for benchmarking.  A negative number will use all reads.
samplerate=<1.0>	Set to a fraction of 1 if you want to randomly sample reads.  For example, samplerate=0.25 would randomly use a quarter of the reads and ignore the rest.  Useful for huge datasets where all you want to know is the % mapped.
sampleseed=<1>		Set to the RNG seed for random sampling.  If this is set to a negative number, a random seed is used; for positive numbers, the number itself is the seed.  Since the default is 1, this is deterministic unless you explicitly change it to a negative number.	
idmodulo=<1>		Set to a higher number if you want to map only every Nth read (for sampling huge datasets).


Mapping Parameters:
fast=<f>		The fast flag is a macro.  It will set many other paramters so that BBMap will run much faster, at slightly reduced sensitivity for most applications.  Not recommended for RNAseq, cross-species alignment, or other situations where long deletions or low identity matches are expected.
minratio=<0.56>		Alignment sensitivity as a fraction of a read's max possible mapping score.  Lower is slower and more sensitive but gives more false positives.  Ranges from 0 (very bad alignment) to 1 (perfect alignment only).  Default varies between BBMap versions. 
minidentity=<>		Or "minid".  Use this flag to set minratio more easily.  If you set minid=0.9, for example, minratio will be set to a value that will be APPROXIMATELY equivalent to 90% identity alignments.
minapproxhits=<1>	Controls minimum number of seed hits to examine a site.  Higher is less accurate but faster (on large genomes).  2 is maybe 2.5x as fast and 3 is maybe 5x as fast on a genome with of gigabases.  Does not speed up genomes under 100MB or so very much.
padding=<4>		Sets extra padding for slow-aligning.  Higher numbers are more accurate for indels near the tips of reads, but slower.
tipsearch=<100>		Controls how far to look for possible deletions near tips of reads by brute force.  tipsearch=0 disables this function.  Higher is more accurate.
maxindel=<16000>	Sets the maximum size of indels allowed during the quick mapping phase.  Set higher (~100,000) for RNA-SEQ and lower (~20) for large assemblies with mostly very short contigs.  Lower is faster.
strictmaxindel=<f>	Set to true to disallow mappings with indels longer than maxindel.  Alternately, for an integer X, 'strictmaxindel=X' is equivalent to the pair of flags 'strictmaxindel=t maxindel=X'.
pairlen=<32000>  	Maximum distance between mates allowed for pairing.
requirecorrectstrand=<t>	Or "rcs".  Requires correct strand orientation when pairing reads.  Please set this to false for long mate pair libraries!
samestrandpairs=<f>	Or "ssp".  Defines correct strand orientation when pairing reads.  Default is false, meaning opposite strands, as in Illumina fragment libraries.  "ssp=true" mode is not fully tested.
killbadpairs=<f>	Or "kbp".  When true, if a read pair is mapped with an inappropriate insert size or orientation, the read with the lower mapping quality is marked unmapped.
rcompmate=<f>		***TODO*** Set to true if you wish the mate of paired reads to be reverse-complemented prior to mapping (to allow better pairing of same-strand pair libraries).
kfilter=<-1>		If set to a positive number X, all potential mapping locatiosn that do not have X contiguous perfect matches with the read will be ignored.  So, reads that map with "kfilter=51" are assured to have at least 51 contiguous bases that match the reference.  Useful for mapping to assemblies generated by a De Bruijn graph assembly that used a kmer length of X, so that you know which reads were actually used in the assembly.
threads=<?>		Or "t".  Set number of threads.  Default is # of logical cores.  The total number of active threads will be higher than this, because input and output are in seperate threads.
perfectmode=<f>		Only accept perfect mappings.  Everything goes much faster.  
semiperfectmode=<f>	Only accept perfect or "semiperfect" mappings.  Semiperfect means there are no mismatches of defined bases, but up to half of the reference is 'N' (to allow mapping to the edge of a contig).
rescue=<t>		Controls whether paired may be rescued by searching near the mapping location of a mate.  Increases accuracy, with usually a minor speed penalty.
expectedsites=<1>	For BBMapPacBioSkimmer only, sets the expected number of correct mapping sites in the target reference.  Useful if you are mapping reads to other reads with some known coverage depth.
msa=<>			Advanced option, not recommended.  Set classname of MSA to use.
bandwidth=0		Or "bw".  When above zero, restricts alignment band to this width.  Runs faster, but with reduced accuracy for reads with many or long indels.
bandwidthratio=0	Or "bwr".  When above zero, restricts alignment band to this fraction of a read's length.  Runs faster, but with reduced accuracy for reads with many or long indels.
usequality=<t>		Or "uq".  Set to false to ignore quality values when mapping.  This will allow very low quality reads to be attempted to be mapped rather than discarded.
keepbadkeys=<f>		Or "kbk".  With kbk=false (default), read keys (kmers) have their probability of being incorrect evaluated from quality scores, and keys with a 94%+ chance of being wrong are discarded.  This increases both speed and accuracy.
usejni=<f>		Or "jni".  Do alignments in C code, which is faster.  Requires first compiling the C code; details are in /jni/README.txt.  This will produce identical output.
maxsites2=<800>		Don't analyze (or print) more than this many alignments per read.
minaveragequality=<0>	(maq) Discard reads with average quality below this.

Post-Filtering Parameters:

idfilter=0              Different than "minid".  No alignments will be allowed with an identity score lower than this value.  This filter occurs at the very end and is unrelated to minratio, and has no effect on speed unless set to 1.  Range is 0-1.
subfilter=-1            Ban alignments with more than this many substitutions.
insfilter=-1            Ban alignments with more than this many insertions.
delfilter=-1            Ban alignments with more than this many deletions.
indelfilter=-1          Ban alignments with more than this many indels.
editfilter=-1           Ban alignments with more than this many edits.
inslenfilter=-1         Ban alignments with an insertion longer than this.
dellenfilter=-1         Ban alignments with a deletion longer than this.

Output Parameters:
out=<outfile.sam>	Write output to this file.  If out=null, output is suppressed.  If you want to output paired reads to paired files, use a "#" symbol, like out=mapped#.sam.  Then reads1 will go to mapped1.sam and reads2 will go to mapped2.sam. (NOTE: split output currently diabled for .sam format, but allowed for native .txt format).  To print to standard out, use "out=stdout"
outm=<>			Write only mapped reads to this file (excluding blacklisted reads, if any).
outu=<>			Write only unmapped reads to this file.
outb=<>			Write only blacklisted reads to this file.  If a pair has one end mapped to a non-blacklisted scaffold, it will NOT go to this file. (see: blacklist)
out2=<>			If you set out2, outu2, outm2, or outb2, the second read in each pair will go to this file.  Not currently allowed for SAM format, but OK for others (such as fasta, fastq, bread).
overwrite=<f>		Or "ow".  Overwrite output file if it exists, instead of aborting.
append=<f>		Or "app".  Append to output file if it exists, instead of aborting.
ambiguous=<best>	Or "ambig". Sets how to handle ambiguous reads.  "first" or "best" uses the first encountered best site (fastest).  "all" returns all best sites.  "random" selects a random site from all of the best sites (does not yet work with paired-ends).  "toss" discards all sites and considers the read unmapped (same as discardambiguous=true).  Note that for all options (aside from toss) ambiguous reads in SAM format will have the extra field "XT:A:R" while unambiguous reads will have "XT:A:U".
ambiguous2=<best>	(for BBSplit only) Or "ambig2". Only for splitter mode.  Ambiguous2 strictly refers to any read that maps to more than one reference set, regardless of whether it has multiple mappings within a reference set.  This may be set to "best" (aka "first"), in which case the read will be written only to the first reference to which it has a best mapping; "all", in which case a read will be written to outputs for all references to which it maps; "toss", in which case it will be considered unmapped; or "split", in which case it will be written to a special output file with the prefix "AMBIGUOUS_" (one per reference).
outputunmapped=<t>	Outputs unmapped reads to primary output stream (otherwise they are dropped).
outputblacklisted=<t>	Outputs blacklisted reads to primary output stream (otherwise they are dropped).
ordered=<f>		Set to true if you want reads to be output in the same order they were input.  This takes more memory, and can be slower, due to buffering in multithreaded execution.  Not needed for singlethreaded execution.
ziplevel=<2>		Sets output compression level, from 1 (fast) to 9 (slow).  I/O is multithreaded, and thus faster when writing paired reads to two files rather than one interleaved file.
nodisk=<f>		"true" will not write the index to disk, and may load slightly faster.   Prevents collisions between multiple bbmap instances writing indexes to the same location at the same time.
usegzip=<f>		If gzip is installed, output file compression is done with a gzip subprocess instead of with Java's native deflate method.  Can be faster when set to true.  The output file must end in a compressed file extension for this to have effect.
usegunzip=<f>		If gzip is installed, input file decompression is done with a gzip subprocess instead of with Java's native inflate method.  Can be faster when set to true.
pigz=<f>          	Spawn a pigz (parallel gzip) process for faster compression than Java or gzip.  Requires pigz to be installed.
unpigz=<f>        	Spawn a pigz process for faster decompression than Java or gzip.  Requires pigz to be installed.
bamscript=<filename>	(bs for short) Writes a shell script to <filename> with the command line to translate the sam output of BBMap into a sorted bam file, assuming you have samtools in your path.
maxsites=<5>		Sets maximum alignments to print per read, if secondary alignments are allowed.  Currently secondary alignments may lack cigar strings.
secondary=<f>		Print secondary alignments.
sssr=<0.95>  		(secondarysitescoreratio) Print only secondary alignments with score of at least this fraction of primary.
ssao=<f>     		(secondarysiteasambiguousonly) Only print secondary alignments for ambiguously-mapped reads.
quickmatch=<f>		Generate cigar strings during the initial alignment (before the best site is known).  Currently, this must be enabled to generate cigar strings for secondary alignments.  It increases overall speed but may in some very rare cases yield inferior alignments due to less padding.
local=<f>          	Output local alignments instead of global alignments.  The mapping will still be based on the best global alignment, but the mapping score, cigar string, and mapping coordinate will reflect a local alignment (using the same affine matrix as the global alignment).
sortscaffolds=<f> 	Sort scaffolds alphabetically in SAM headers to allow easier comparisons with Tophat (in cuffdif, etc).  Default is in same order as source fasta.
trimreaddescriptions=<f>	(trd) Truncate read names at the first whitespace, assuming that the remaineder is a comment or description.
machineout=<f>    	Set to true to output statistics in machine-friendly 'key=value' format.
forcesectionname=<f>	All fasta reads get an _# at the end of their name.  The number is 1 for the first shred and continues ascending.


Sam settings and flags:
samversion=<1.4>  	SAM specification version. Set to 1.3 for cigar strings with 'M' or 1.4 for cigar strings with '=' and 'X'.  Samtools 0.1.18 and earlier are incompatible with sam format version 1.4 and greater.
saa=<t>           	(secondaryalignmentasterisks) Use asterisks instead of bases for sam secondary alignments.
cigar=<t>		Generate cigar strings (for bread format, this means match strings).  cigar=false is faster.  "cigar=" is synonymous with "match=".  This must be enabled if match/insertion/deletion/substitution statistics are desired, but the program will run faster with cigar strings disabled.
keepnames=<f>		Retain original names of paired reads, rather than ensuring both reads have the same name when written in sam format by renaming read2 to the same as read1.  If this is set to true then the output may not be sam compliant.
mdtag=<f>		Generate MD tags for SAM files.  Requires that cigar=true.  I do not recommend generating MD tags for RNASEQ or other data where long deletions are expected because they will be incredibly long.
xstag=<f>		Generate XS (strand) tags for Cufflinks.  This should be used with a stranded RNA-seq protocol.
xmtag=<t>		Generate XM tag.  Indicates number of best alignments.  May only work correctly with ambig=all.
nhtag=<f>             	Write NH tags.
intronlen=<999999999>	Set to a lower number like 10 to change 'D' to 'N' in cigar strings for deletions of at least that length.  This is used by Cufflinks; 'N' implies an intron while 'D' implies a deletion, but they are otherwise identical.
stoptag=<f>		Allows generation of custom SAM tag YS:i:<read stop location>
idtag=<f>		Allows generation of custom SAM tag YI:f:<percent identity>
scoretag=<f>		Allows generation of custom SAM tag YR:i:<raw mapping score>
inserttag=<f>		Write a tag indicating insert size, prefixed by X8:Z:
rgid=<>     		Set readgroup ID.  All other readgroup fields can be set similarly, with the flag rgXX=value.
noheader=<f>		Suppress generation of output header lines.


Statistics and Histogram Parameters:
showprogress=<f>	Set to true to print out a '.' once per million reads processed.  You can also change the interval with e.g. showprogress=20000.
qhist=<file>		Output a per-base average quality histogram to <file>.
aqhist=<file>		Write histogram of average read quality to <file>.
bqhist=<file>		Write a quality histogram designed for box plots to <file>.
obqhist=<file>		Write histogram of overall base counts per quality score to <file>.
qahist=<file>		Quality accuracy histogram; correlates claimed phred quality score with observed quality based on substitution, insertion, and deletion rates.
mhist=<file>		Output a per-base match histogram to <file>.  Requires cigar strings to be enabled.  The columns give fraction of bases at each position having each match string operation: match, substitution, deletion, insertion, N, or other.
ihist=<file>		Output a per-read-pair insert size histogram to <file>.
bhist=<file>		Output a per-base composition histogram to <file>.
indelhist=<file>  	Output an indel length histogram.
lhist=<file>      	Output a read length histogram.
ehist=<file>      	Output an errors-per-read histogram.
gchist=<file>     	Output a gc content histogram.
gchistbins=<100>	(gcbins) Set the number of bins in the gc content histogram.
idhist=<file>    	Write a percent identity histogram.
idhistbins=<100>	(idbins) Set the number of bins in the identity histogram.
scafstats=<file>	Track mapping statistics per scaffold, and output to <file>.
refstats=<file>		For BBSplitter, enable or disable tracking of read mapping statistics on a per-reference-set basis, and output to <file>.
verbosestats=<0>	From 0-3; higher numbers will print more information about internal program counters.
printunmappedcount=<f>	Set true to print the count of reads that were unmapped.  For paired reads this only includes reads whose mate was also unmapped.


Coverage output parameters (these may reduce speed and use more RAM):
covstats=<file>		Per-scaffold coverage info.
covhist=<file>		Histogram of # occurrences of each depth level.
basecov=<file>		Coverage per base location.
bincov=<file>		Print binned coverage per location (one line per X bases).
covbinsize=1000		Set the binsize for binned coverage output.
nzo=f			Only print scaffolds with nonzero coverage.
twocolumn=f		Change to true to print only ID and Avg_fold instead of all 6 columns to the 'out=' file.
32bit=f			Set to true if you need per-base coverage over 64k.
bitset=f		Store coverage data in BitSets.
arrays=t		Store coverage data in Arrays.
ksb=t			Keep residual bins shorter than binsize.
strandedcov=f		Track coverage for plus and minus strand independently.  Requires a # symbol in coverage output filenames which will be replaced by 1 for plus strand and 2 for minus strand.
startcov=f		Only track start positions of reads.
concisecov=f		Write basecov in a more concise format.


Trimming Parameters:
qtrim=<f>		Options are false, left, right, or both.  Allows quality-trimming of read ends before mapping.
			false: Disable trimming.
			left (l): Trim left (leading) end only.
			right (r): Trim right (trailing) end only.  This is the end with lower quality many platforms.
			both (lr): Trim both ends.
trimq=<5>		Set the quality cutoff.  Bases will be trimmed until there are 2 consecutive bases with quality GREATER than this value; default is 5.  If the read is from fasta and has no quality socres, Ns will be trimmed instead, as long as this is set to at least 1.
untrim=<f>		Untrim the read after mapping, restoring the trimmed bases.  The mapping position will be corrected (if necessary) and the restored bases will be classified as soft-clipped in the cigar string.


Java Parameters:
-Xmx       		If running from the shellscript, include it with the rest of the arguments and it will be passed to Java to set memory usage, overriding the shellscript's automatic memory detection.  -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max allowed is typically 85% of physical memory.
-da			Disable assertions.  Alternative is -ea which is the default.


Splitting Parameters:
The splitter is invoked by calling bbsplit.sh (or align2.BBSplitter) instead of bbmap.sh, for the indexing phase.  It allows combining multiple references and outputting reads to different files depending on which one they mapped to best.  The order in which references are specified is important in cases of ambiguous mappings; when a read has 2 identically-scoring mapping locations from different references, it will be mapped to the first reference.
All parameters are the same as BBMap with the exception of the ones listed below.  You can still use "outu=" to capture unmapped reads.
ref_<name>=<fasta files>	Defines a named set of organisms with a single fasta file or list.  For example, ref_a=foo.fa,bar.fa defines the references for named set "a"; any read that maps to foo.fasta or bar.fasta will be considered a member of set a.
out_<name>=<output file>	Sets the output file name for reads mapping to set <name>.  out_a=stuff.sam would capture all the reads mapping to ref_a.
basename=<example%.sam>		This shorthand for mass-specifying all output files, where the % symbol is a wildcard for the set name.  For example, "ref_a=a.fa ref_b=b.fa basename=mapped_%.sam" would expand to "ref_a=a.fa ref_b=b.fa out_a=mapped_a.sam out_b=mapped_b.sam"
ref=<fasta files>		When run through the splitter, this is shorthand for writing a bunch of ref_<name> entries.  "ref=a.fa,b.fa" would expand to "ref_a=a.fa ref_b=b.fa".


Formats and Extensions
.gz,.gzip,.zip,.bz2	These file extensions are allowed on input and output files and will force reading/writing compressed data.
.fa,.fasta,.txt,.fq,.fastq	These file extensions are allowed on input and output files.  Having one is REQUIRED.  So, reads.fq and reads.fq.zip are valid, but reads.zip is NOT valid.  Note that outputting in fasta or fastq will not retain mapping locations.
.sam			This is only allowed on output files.
.bam			This is allowed on output files if samtools is installed.  Beware of memory usage; samtools will run in a subprocess, and it can consume over 1kb per scaffold of the reference genome.


Different versions:
BBMap			(bbmap.sh) Fastest version.  Finds single best mapping location.
BBMapPacBio		(mapPacBio.sh) Optimized for PacBio's error profile (more indels, fewer substitutions).  Finds single best mapping location.  PacBio reads should be in fasta format.
BBMapPacBioSkimmer	(bbmapskimmer.sh) Designed to find ALL mapping locations with alignment score above a certain threshold; also optimized for Pac Bio reads.
BBSplitter		(bbsplit.sh) Uses BBMap or BBMapPacBio to map to multiple references simultaneously, and output the reads to the file corresponding to the best-matching reference.  Designed to split metagenomes or contaminated datasets prior to assembly.
BBWrap			(bbwrap.sh) Maps multiple read files to the same reference, producing one sam file per input file.  The advantage is that the reference/index only needs to be read once.




Notes.
File types are autodetected by parsing the filename.  So you can name files, say, out.fq.gz or out.fastq.gz or reads1.fasta.bz2 or data.sam and it will work as long as the extensions are correct.


Change Log:

V34.
34.00
Fixed a bug in BandedAlignerConcrete related to width being allowed to be even.
34.01
IdentityMatrix is now much faster for ghigh-identity sequences, and allows the 'width' flag to increase speed.
Updated FilterReadsByName to allow "names=<read filename>", supporting fastq, fasta, and sam.  So, one file will be filtered according to the names of reads in a second file.  "names=<file>" where the file is just a list of names is still supported.
34.02
Fixed a couple errors in ConcurrentReadInputStreamD.
Added fetching of a dummy list from "empty" for crisD, both master and slave.
Added A_SampleD, which uses crisD.  It now works correctly for master.
Renamed various ConcurrentReadStreamInterface classes.
Added an abstract superclass for all ConcurrentReadInputStreams, which extends Thread.  Now, cris can be started directly without making a new thread.
Changed all instances of wrapping cirs in a thread to just use start directly.  These are mostly commented with "//4567" to find if something was missed (like starting the cris twice).
Increased cris stability by removing "returnList(ListNum, boolean)" and replacing it with "returnList(long, boolean)".  Lists may no longer be recycled.
34.03
Added scaffoldstats to BBQC and RQCFilter fileList logs.  Requested by Bryce F.
Fixed a strange deadlock in Dedupe/ConcurrentCollectionReadInputStream caused by making CRIS a Thread subclass.  This will still occur if CRIS goes back to being a Thread.  Noted by Shoudan.
34.04
Removed hitCount tracking from Seal.
"qtrim=<integer>" is now allowed for all classes using Parser.parseTrim().
Parser.parseZip, parseInterleaved, parseQuality, parseTrim, parseFasta, and parseCommonStatic were integrated into most classes; reduced code size by almost 200kb.
Parser.parseTrim got some extra functionality, like maxNs.
Made an abstract superclass for KmerCount* classes, allowing removal of some code.
Removed all KmerCount.countFasta() methods; they must now use a CRIS.
Retired ErrorCorrectMT (superceded ny KmerNormalize).
Fixed bug in BBDuk, Seal, and ReformatReads - when quality trimming and force-trimming, count of trimmed reads could go over 100%.  Now these counts are independent.  Noted by ysnapus (SeqAnswers).
Removed "minscaf" and "mincontig" flags from Parser.parseFasta() because they were conflated.
Determined cause of Kurt's error message in Dedupe - lower-case letters can trigger a failure.
Dedupe now defaults to "tuc=t" (all input is made upper-case).
Moved CRIS factory from CGRIS to CRIS.
Copied cc2-cc5 to /global/projectb/sandbox/gaag/TestData/SingleCell/SimMockCommunity/plate*/.  These are simulated cross-contaminated single cell plates.
Removed conflated "qual" flag from RandomReads; "q" should be used instead to set all read quality values to a single number.
Fixed conflated "renamebymapping" flag in RenameReads.
"tbr" flag is conflated in KmerNormalize; adjusted so that it now controls both "tossBadReads" (reads with errors) and "tossBrokenReads" (reads with the wrong number of quality scores).
Conflated "gzip" flag in ChromArrayMaker/FastaToChromArrays changed to "gz".
Handled conflated "ziplevel" flag in AbstractMapper.
Conflated "fakequality" flag resolved by moving from BBMap to Parser and renaming "fakefastaquality"/"ffq".
Added hdist2 and edist2 to BBDuk.  These allow independently specifying hdist/edist for full-length kmers and short kmers when using mink.
Added trimhdist2 to RQCFilter/BBQC.
*** Added path and mapref flags to RQCFilter/BBQC; they can now map to an arbitrary genome instead of just human.
Added Shared.USE_MPI field (parsed by Parser.parseCommonStatic; "mpi" or "usempi").
Added Shared.MPI_RANK field (should be set automatically).
Added Shared.MPI_KEEP_ALL field.  This controls whether CRISD objects retain all reads, or just some of them.
CRIS now automatically returns a CRISD when USE_MPI is enabled, as a slave or master depending on whether rank==0.
ListNum is now Serializable.
CRISD now transmits ListNum objects rather than ArrayLists, so that the number is preserved.
Added Maxns flag to reformat.
Fixed BBQC and RQCFilter's unnecessary addition of "usejni" to BBMap phase, since it is now already parsed by parseCommonStatic.
BBQC now defaults to normalization and ecc off, but can be enabled with the "ecc" and "norm" flags, and supports cecc flag.
Added notes on compiling JNI version suggested by sdriscoll.
34.05
Commented out a reference to ErrorCorrectMT in MateReadsMT.
34.06
FindPrimers (msa.sh) now accepts multiple queries (primers) and will use the best-matching of them.
Added a BBMap flag to disable score penalty due to ambiguous alignments ("penalizeambiguous" or "pambig").  Requested by Matthias.
Fixed failure to start CRIS in A_SampleD.
Fixed some incorrect division in CRISD.
Added MPI_NUM_RANKS to Shared.  This is parsed by parser via e.g. "mpi=4".
Added BBMap flags subfilter, insfilter, delfilter, inslenfilter, dellenfilter, indelfilter, editfilter.  These function similarly to idfilter.  Requested by sdriscoll.
34.07
Dedupe now automatically calls Dedupe2 if more than 2 affixes are requested.
Added "subset" (sst) and "subsetcount" (sstc) flags to Dedupe.
Added "printLengthInEdges" (ple) flag to Dedupe.
34.08
Finished Dedupe subset processing for graph file generation.
34.09
Fixed bug where 'k' was not added to filename in RQCFilter.  Noted by Vasanth.
34.10
Documented "ordered" and "trd" flags for BBDuk/Seal.
Added crismpi flag to allow switching between crisd and crismpi.
Added shared.mpi package, containing MPIWrapper and ConcurrentReadInputStreamMPI.
34.11
Added detection of read length, interleaving, and quality coding to FileFormat objects, but these fields are not currently read.
FileFormat.main() now outputs read length, if in fastq format.
Reformat will now allow sam -> sam conversion; not useful in practice, but maybe useful in testing.
Added flag "mpikeepall", default true.
Fixed deadlock when mpikeepall=false.  Noted by Jon Rood.
34.12
Added 'auto' option to gcbins and idbins flags.  Requested by Seung-jin.
Added dedupe "addpairnum" flag to control whether ".1" is appended to numbered graph nodes.
Added real quality to qhist plot, when mhist is being generated.
Moved maxns and maq to AFTER quality trimming in RQCFilter and BBDuk.
Added "ftm" (forcetrimmodulo) flag to BBDuk/Reformat/RQCFilter/BBQC.  Default 5 for RQCFilter/BBQC, 0 otherwise.
34.13
Fixed a missing "else" in RQCFilter/BBQC.  Noted by Kurt LaButti.
34.14
Added .size() to ListNum.
CrisD gained "unicast" method.  Also, unicast and listen now have mandatory toRank parameter.
Made CrisD MPI methods protected rather than private, so they can be overridden.
Refactored RTextOutputStream3.
34.15
Added Shared.LOW_MEMORY:
Disables multithreaded index gen.
Disables multithreaded ReadWrite writeObjectInThread method.
Disables ByteFile2.
For some reason it does not really seem to reduce memory consumption...
Added BBMap qfin1 and qfin2 flags.
Updated BBMap to use more modern input stream initialization.
Added mapnt.sh for mapping to nt on a 120g node.
34.16
Changed RQCFilter "t" to mean "trimmed"; "k" was removed.
Added parser noheadersequences (nhs) flag for sam files with millions of ref sequences.
Documented "ambig" flag in Seal.
Fixed issue where Shared.READ_BUFFER_NUM_BUFFERS was not getting changed with THREADS was changed.  Now both are private and get set together.
Verified that mapnt.sh works on 120G nodes.
34.17
RTextOutputStream3 renamed to ConcurrentReadOutputStream.
ReadStreamByteWriter refactored to be cleaner.
Merged MPI dev branch into master.
34.18
Moved Seal's maxns/maq to after trimming.
Added chastity filter to bbduk and reformat (reads containing " 1:Y:" or " 2:Y:").  Requested by Lynn A.
Dedupe outd stream now produces correctly interleaved reads.  Requested by Lynn A.
Replaced Dedupe TextStreamWriters with ByteStreamWriters, for read output.
34.19
Added parseCommon() to BBDuk, allowing samplerate flag.
34.20
FASTA_WRAP moved to Shared.
Numeric qual output is now wrapped after the same number of bases as fasta output.
"Low quality discards:" line is now triggered by chastity filter.
SPLIT_READS and MIN_READ_LEN are now disabled when processing reference in BBDuk/Seal.
Seal gained parseCommon and parseQuality.
34.21
Fixed MIN_READ_LEN bug (set to 0; should have been set to 1)
34.22
Added qfin (qual file) flags to BBDuk/Seal.
Applied BBDuk restrictleft and restrictright to filtering and masking; before, it was only valid for trimming.
Added calcCigarBases.
Required includeHardClip parameter for all calls to calcCigarLength(), start(), or stop().
Fixed bug in pileup caused by hard-clipped reads.  Noted by Casey B.
34.23
DecontaminateByNormalization was excluding contigs with length under 50bp, which caused an assertion error.
Fixed a crash in BBDuk2 when not using a reference.  Noted by Dave O.
Added entropy filter to BBDuk/BBDuk2.  Set "entropy=X" where X is 0 to 1 to filter reads with less entropy than X.
34.24
Added maxreads flag to readlength.sh.
Fixed bug in BBMap - when directly outputting coverage, secondary alignments were never being used.
BBMap now uses the "ambig" and "secondary" flags to determine whether to include secondary site coverage.  Specifically, "ambig=all" will use secondary sites, while other modes will not unless "secondary=t".  In other words, use of secondary sites in coverage will be exactly the same as use of them in a sam output file.  Removed "uscov=t                 Include secondary alignments when calculating coverage." from shellscript.
Fixed minid trumping minratio when both were specified.  Now, the last one specified will be used.
Added pileup support for reads with asterisks instead of bases, as long as they have a cigar string.  Also sped up calculation of read stop position.
Cigar string 'M' symbols are now converted to match string 'N' symbols if there is no reference.
34.25
BBMerge initialization order bug fixed; it was preventing jni from being used with the "loose" or "vloose" flags.  Noted by sarvidsson (SeqAnswers).
34.26
Fixed semiperfect mode allowing non-semiperfect rescued alignments.  Noted by Dave O.
Fixed ReadStats columns header for qhist when mhist was also generated.
Fixed an inequality in BBMergeOverlapper that favored shorter overlaps with an equal number of mismatches, in some cases.  Had no impact on a normal 1M read benchmark except when margin=0, where it tripled the false-positive rate.
34.27
Enabled verbose mode in BBMergeOverlapper.
34.28
Added "align2." to sam header command line of BBMap.
Fixed bug in BBMap that could cause "=" to be printed for "rnext" even when pairs were on different scaffolds.  Noted by rkanwar (SeqAnswers).
34.30
Reformat can now produce indelhist from sam files prior to v1.4.
Fixed a crash bug in BBMap caused by an improper assertion.  Noted by Rob Egan.
34.31
BBDuk/Seal now recognize "scafstats" flag as equivalent to "stats".
Seal now defaults to 5 stats columns (includes #bp).
Wrote BBTool_ST, and abstract superclass for singlethreaded BBTools.
Clarified documentation of "trimq=X" as meaning "regions with average quality under X will be trimmed".
Fixed major bug in RQCFilter/BBQC:  "forcetrimmod" was being set to same value as "ktrim".  Noted by Brian Foster.
34.32
Changed the way BBMerge handles qualities to make it 40% faster (in java mode).  Reduced size of jni matrix accordingly.
Fixed lack of readgroup tags for unmapped reads in sam format.  Noted by Rahul (SeqAnswers).
Ensured Read.CHANGE_QUALITY affects both lower (<0) and upper (>41) values.
34.33
Pushed BBMergeOverlapper.c to commit.
34.34
Documented trimfragadapter and removehuman in RQCFilter.
Added Parser flag for Shared.READ_BUFFER_LENGTH (readbufferlength).
Added Parser flag for Shared.READ_BUFFER_MAX_DATA (readbufferdata).
Added Parser flag for Shared.READ_BUFFER_NUM_BUFFERS (readbuffers).
RQCFilter now accepts multiple references for decontamination by mapping.
Added FuseSequence (the first BBTool_ST subclass) and fuse.sh, for gluing contigs together with Ns.
Reformatted many scripts' help info to remove echo statements.
Fixed bugs in stats and countgc; they were not including undefined bases when printing the length in gcformat=1 and gcformat=4.
Replaced all instances of .bases.length with .length(), to prevent null pointer exceptions (for example in sam lines with no bases).
Added cat and dog flags to rqcfilter.
Changed defaults of BBMask to reduce amount masked in cat and dog to ~1% of genome.  This still masks all of the coincidental low-complexity hits from fungi.
Determined that dog is contaminated with fungus, particularly chr7 and chr13.
34.35
Fixed a bug in which data was retained from the prior index when indexing a second fasta file in nodisk mode.
34.36
Disabled an assertion in BBMerge that the input is paired; it crashes if the input file is empty.
34.37
NSLOTS is now ignored if at least 16, to account for new 20-core nodes. 
ReadWrite.getOutputStream now creates the directory structure if it does not already exist.  Problem discovered by Brian Foster.
BBQC and RQCFilter now strip directory names before writing temp files.
BBDuk now correctly reports number of reads quality-filtered.
Added "unmappedonly" flag to Reformat.
RQCFilter now defaults to using TMPDIR.
34.38
BBMap now prints reads/second correctly.  Before, it actually displayed pairs/second with paired data.
Added maxq flag to BBMerge, which allows quality values over 41 where reads overlap.  Requested by Eric J.
Changed CoveragePileup from TextFiles to ByteFiles; increased read speed by 3.65x.
Changed CoveragePileup from TextStreamWriters to ByteStreamWriters; increased write speed by 1.46x.
Fixed a bug in BBQC/RQCFilter: paired input and interleaved output was getting its paired status lost.  Noted by Simon P.
Reformat, when in "mappedonly" or "unmappedonly" mode, now excludes reads with no bases or secondary alignments.


TODO: BED support for pileup.  And make Pileup faster by ignoring irrelevent sam fields.
TODO: CrossMask.  Accept set of files; for each, mask using BBDuk with all others as ref.
TODO: Study bisulfite data on BBMap.  Possibly use multiple reference copies with different transforms (C->T, A-G, both, neither).  
TODO: Shellscripts are not able to handle paths containing spaces.
TODO: Add mininsert flag to BBMap.  And maybe maxinsert.
TODO: Send S.L. CC reads only.
TODO: Parse MD tag when available.
TODO: CC rates for all 3 platforms in one chart; ignore R1/R2 differences.
TODO: Add BBDuk/Seal option to mutate read kmers, not ref kmers.
*TODO: Dedupe loses reads when using paired data and run multithreaded.
*TODO: strange bug in BBMap causing reads to be mapped to different chromosomes but display "=" for ref sequence.  Can replicate with GATK human ref.
TODO: document nhs flag.
TODO: Filter cross-contam plates with only depth and length, test cc rates.
TODO: Fix dedupe crash when minclustersize=1.
TODO: Clarify or fix what minid does in Dedupe.
TODO: Finish cat+dog.
TODO: Add ribosomal filtering to rqc.
TODO: Update BandedAlignerJNI for quicker width reset.
TODO: Optional penalty when seq ends before ref in banded.
TODO: Make sure AddAdapters is adding them correctly, i.e., reverse-complemented (or not).
TODO: Make list of proposed higher stringency adapter trimming changes and send to Vasanth/Erika.
TODO: Retire ErrorCorrect, and move the functionality over to another class.
TODO: Implement ErrorCorrectBulk in KmerNormalize.  It is used in MateReadsMT.
TODO: BBMerge should allow optional inline error-correction for reads that fail to merge, and revert if they still fail.
TODO: Retire KmerCount7MT (non-atomic version).
*TODO: It appears that timeslip is being correctly applied by fillLimited (etc), but not by calcDelScore() or calcAffineScore().
TODO: Dedupe should warn if lowercase letters are present. (Kurt)


v33.
Added "usemodulo" flag to BBMap.  Allows throwing away 80% of reference kmers to save memory.  Slight reduction in sensitivity.  Requested by Rob Egan.
Moved GetReads back to jgi package and fixed shellscript.
Fixed rare crash when using "local" mode on paired-end data on highly-repetitive genomes (Creinhardtii).  Found by Vasanth S.
Improved "usemodulo" mode - it was biased against minus-strand hits.  Now, it keeps kmers where (kmer%5==rkmer%5).  Result is virtually no reduction in sensitivity (zero in error-free reads, and less than 0.01% in reads with 8% error).
BBMap will now discard reads shorter than "minlen".
Added "idhistbins" or "idbins" flag to BBMap; allows setting the number of bins used in the idhist.
Rescaled BBMap's MAPQ to be lower.  It is now 0 for unmapped, 1-3 for ambiguous, and roughly 4-45 otherwise, with higher values allowed for longer reads.
Added a much flatter MSA version, "MultiStateAligner9Flat", requested by JJ Chai.
Fixed SNR output formatting.
Added "forcesectionname" flag; fasta reads will always get an "_1" at the end, even if they are not broken into multiple pieces. (requested by Shoudan)
Changed "fastareadlen" suffixes to only be appended when read is > maxlen rather than >=
Reorganized SamLine and created SamHeader class.
Modified CountBarcodes to append sub distance from expected barcodes and 'valid' for valid barcodes.
Fixed null pointer exception related to "qhist", "aqhist", and "qahist".  Noted by Harald (seqanswers).
Fixed issue of readlength.sh breaking up reads when processing fasta files without a fasta extension.
Updated BBDuk documentation.
Added "maxlength" and qahist support to BBDuk.
Added "minoverlap" and "mininsert" to BBDuk.
Added "maxlength" to BBMerge.
Created countbarcodes.sh
Added edit distance column to CountBarcodes output.
Added raw mapping score tag, YS:i:, controlled by "scoretag" flag and disabled by default.
Added 'cq' (changequality) flag to reformat.  Default: true.
Fixed mhist being generated from sam files.
Added readgroup support; a readgroup field "xx" can be specified with the flag "rgxx=value".
Updated 'usemodulo' flag to use (kmer%9==0 || rkmer%9==0).  Requiring the remainders to be equal unevenly affected palindromes and thus even kmer lengths.
Updated RemoveHuman to use 'usemodulo' flag and reduced RAM allotment from 23g to 10g.  Updated index location of HG19 masked.
Added "idfilter" to BBMap.
Made BandedAligner abstract superclass and created BandedAlignerConcrete for the Java implementation, and BandedAlignerJNI for the C version.
Made file extension detection more robust against capitalization.
Added outsingle to BBDuk.
Replaced FastaToChromArrays with ChromArrayMaker.  Now, indexing can be done from fastq files instead of just fasta.
Fixed MAJOR bug in which reference was split up into pieces (as of 33.12).
Reverted to old version of reference loader (as of 33.13) as there was still a bug (skipping every other scaffold).
BBDuk (and BBDuk2) now better support kmer masking!  Every occurance of a kmer is individually masked.
Added parseQuality (qin, qout, etc) to Dedupe.
Changed Dedupe default cluster stats cutoff to 2 (from 10), min cluster size to 2, and by default these values are linked.
Added 'outbest' to Dedupe, writing the representative read per cluster (regardless of 'pbr' flag).  This is mainly for 16s clustering.
Fixed sorting of depths in pileup.sh.  Noted by Alicia Clum.
Fixed 'outbest' of Dedupe (was writing to wrong stream).
Slightly accelerated read trimming.
Added read/base count tracking to ConcurrentReadStreamInterface.
Added display of exact number of input and output bases and reads to reformat.sh (requested by Seung-Jin).
Fixed capital letters changing to lower-case in output filenames when using the "basename" flag with BBSplit.  Noted by Shoudan Liang.
Added Tools.condenseStrict(array).
Fixed fast/slow flags with BBSplit.  Noted by Shoudan Liang.
Added 3-frames option to TranslateSixFrames by adding the flag "frames=3".  Requested by Anne M.
TranslateSixFrames now defaults to fasta format when the file extension is unclear.
Added "estherfilter.sh" for filtering blastall queries.
Added option of getting an input stream from a process with null file argument.
Wrote FastaToChromArrays2 based on ByteFile/ByteBuilder for slightly better indexing speed and lower memory use.
Modified ChromosomeArray to work with ByteBuilder.
Fixed reformat displaying wrong number of input reads when run interleaved (due to recent changes).
Added minratio, maxindel, minhits, and fast flags to BBQC, for controlling BBMap.
Fixed "assert(false)" statement accidentally left in SamPileup from testing.  Noted by Brian Foster.
Added kfilter and local flags to BBQC.
Fixed "bs" (bamscript) flag with BBSplit.  Previously, it did not include the per-reference output streams.
Added Jonathan Rood's C code and JNI class for Dedupe.
Modified dedupe shellscripts to allow JNI code.
BBSplit was not outputting any reads when reference files had uppercase letters (as a result of the recent case-sensitivity change).  This has been fixed.  Noted by Shoudan Liang.
BBMap can now output fastq files with reads renamed to indicate mapping location, using the flags "rbm" and "don" (renamebymapping and deleteoldname).
FastaQualInputStream replaced by FastaQualInputStream3.  At least 2.5x faster, and correctly reads input in which fasta and qual lines are wrapped at different lengths.  Bug noted by Kurt LaButti.
Added bqhist, which allows box plots of read quality-per-base-location.
Fixed a slowdown when making quality histograms due to recalculating probability rather than using cached value.
Default sam format is now 1.4.
RemoveHuman/BBQC/RQCFilter now default to minhits=1 because 'usemodulo' reduces the number of valid keys.
Programs no longer default to outputting to stdout when "out=" is not specified because it's annoying.  To write to stdout set "out=stdout.fq" (for example).
AssemblyStats now counts IUPAC and invalid characters seperately.  X and N now denote gaps between contigs, but no other symbols do.  The code was also cleaned somewhat.  The output formatting changed slightly.
Preliminarily integrated Jon Rood's JNI versions of BandedAligner and MultiStateAligner into both Java code and shellscripts to test Genepool deployment.
C code is now in /jni/ folder, at same level as /resources/ and /docs/.
Clarified documentation of BBMap, BBSplit, and BBWrap to differentiate some parameters.  For example, "refstats" only works with BBSplit.
Added LW and RW (whisker values) columns to bqhist output, set at the 2nd and 98th percentiles.  Requested by Seung-Jin Sul.
BBQC will now compress intermediate files to level 2 instead of level 4, to save time.
Fixed incompatibility of dot graph output and other output in Dedupe.
Reverted to default "minhits=2" for RemoveHuman, because minhits=1 took 5x as long.
Added median, mean, and stdev to gchist.  Requested by Seung-Jin.
Added obqhist (overall base quality histogram).  Requested by Seung-Jin.
Fixed various places, such as BBDuk, where the "int=true" flag caused references to be loaded interleaved.  Noted by Jessica Jarett.
Added some parser flags to allow dynamically enabling verbose mode and assertions specifically for certain classes. 
Fixed a bug in BBMap that made secondary alignments sometimes not get cigar strings.
Added "addprefix" mode to rename reads, which simply prepends a prefix to the existing name.
Clarified documentation of different histogram outputs in shellscripts.
Ported BBMapThread changes over to BBMap variants.
Restructured SamPileup and renamed it to CoveragePileup.  Now supports Read objects (instead of just SamLines).
Integrated CoveragePileup with BBMap and documented new flags.
CoveragePileup: Added a concise coverage output, stranded coverage, and read-start-only coverage.
Removed an obsolete Java classes and some shellscripts.
Increased robustness of BBDuk's detection invalid file arguments, and clarified the error messages.  Noted by Scott D.
Fixed a problem with interleaving not being forced on fasta input.
Paired output files will now force BBDuk input to be treated as interleaved.
BBDuk now tracks statistics on which reference sequences were trimmed or masked - previously, it just tracked what was filtered.
Reverse-complemented Nextera adapters and added them to official release (/resources/nextera.fa.gz).
Added Illumina adapter sequence legal disclaimer to /docs/Legal_Illumina.txt
Implemented GC calculation from index, for generating coverage stats while mapping.
Tracked down strangeness with BBDuk.  It is possible for "rcomp=f" to slightly reduce sensitivity when "mm=t" using an even kmer length, due to asymmetry.  This appears to be correct.
Merged in revised JNI Dedupe version that should be working correctly.  Verified that it returns same answer as non-JNI version.  Tests indicate roughly triple speed, when working with PacBio reads of insert.
BBMap JNI version now seems roughly 30% faster than Java version.
Added insert size quartiles to BBMap and BBMerge.  Requested by Alex Copeland.
Fixed rare bug related to SiteScore.fixXY(), caused by aligning reads with insufficient padding, fixing the tips, but not changing the start/stop positions.  Found by Brian Foster.
Fixed a race condition in TextStreamWriter that could randomly cause a deadlock in numerous different programs.  Found by Shoudan Liang.
Added "maxsites2" flag to allow over 800 alignments for a given read.
Fixed bounds of kmer masking in BBDuk; they were off by 2 (too big).
Fixed unintended debug print line.  Noted by Shoudan Liang.
Updated RandomReadInputStream to work with the newer RandomReads3 class.
ConcurrentGenericReadInputStream now supports RandomReadInputStream3 as a producer.
Fixed kmer dumping from CountKmersExact.
Fixed length of vector created in BBMergeOverlapper (4->5).  Noted by Jon Rood.
Changed default kmer length in BBDuk to 27 so that the 'maskmiddle' base will be in the middle for both forward and reverse kmers.
"pairlen" flag accidentally deleted from BBMap; restored.  Noted by HGV (seqanswers).
BBMerge now has a JNI version from Jonathan Rood - 60% faster than pure Java.  Requires compiling the C code; details are in /jni/README.txt.
Wrapped BBMerge JNI initializer in a conditional, so it will not try to load unless "usejni" is specified.  
Added "parseCommonStatic" to BBMerge and BBDuk (to allow JNI flag parsing).
Commented out "module load" and "module unload" statements in public version.
Added 'printlastbin' or 'plb' flag to countunique to produce a final bin smaller than binsize.  Suggested for use in cumulative mode.  Requested by Andrew Tritt.
Added support for bzip2 and pbzip2 compression and decompression.  The programs must be installed to use bz2 format.
Elminated use of "sh" when launching subprocesses.  This also allows pigz compression support in Windows.
Files were not being closed after "testInterleaved()".  Fixed.
Improved error messages when improper quality values are observed.
Updated hard-coded adapter path to include Nextera adapters.  This affects BBQC and RQCFilter.
Improved file format detection.  Now FileFormat (testformat.sh) will print a warning when the contents and extension don't match, and it can differentiate between sam and fastq.  Problem noted by Vasanth Singan.
Fixed issue where "scafstats" output was printing inflated numbers with chimeric paired reads, or pairs with only one mapped read.  Noted by HGV (seqanswers).
Closed stream after reading in FileFormat.
Unrolled, debranched, and removed assertion function calls from BBMerge inner loop.
Fixed a bug in which findTipDeletions was not changing the bounds of the gap array.
Added getters and setters for SiteScores that enforce gap correctness.
Improved GapTools to test for and fix non-ascending points.
Forced use of setters in TranslateColorspaceRead, AbstractMapThread, and BBIndex* classes; this caught some inconsistencies that should increase stability and correctness.
Enabled jni-mode alignment by default for BBQC and removehuman.
Added a BBMap output line indicating how many reads survived for use with, e.g., removehuman.  Requested by Brian Foster.
Added messages to BBQC to indicate which phase is executing.  Requested by Brian Foster.
SiteScore start and stop are exclusively set by methods now.  Fixed a bug with local flag noted by Vasanth Singan.
Added MaximumSpanningTree generation to Dedupe (mst flag).
Merged in faster BBMerge overlapper JNI version; now 90% faster than Java with fastq and 70% faster with fasta.
Improved Dedupe's support for paired reads: fixed an assertion, and added "in1" and "in2".
Fixed a assertion involving semiperfect alignments of repetitive reads, that go out of the alignment window.  Found by Alicia Clum.
Fixed idhist mean calculation.  Added mode, median, stdev, both by read count and base count.
Better documented ConcurrentReadStreamInterface.
Fixed a crash in CoveragePileup when using 32-bit mode.
Fixed a couple instances in which the first two arguments being unrecognized would not be noticed.
Fixed a bug in pileup causing coverage fraction to be reported incorrectly, if arrays were not being used.  Noted by Vasanth Singan.
Fixed a twocolumn mode in pileup; it was generating no output.
Added additional parse flags to pileup, such as "stats" and "outcov".
Added additional output fields to coverage stats - total number of covered bases, and number of reads mapped to plus and minus strands.
CountKmersExact: Added preallocation (faster, less memory) and a one-pass-mode for the prefilter (faster, but nondeterministic).
Replaced most instances of "Long.parseLong" with "Tools.parseKMG" to support kilo, mega, and giga abbreviated suffixes.
Added jgi.PhylipToFasta and phylip2fasta.sh, for converting interleaved phylip files to fasta.  Requested by Esther Singer.
v33.58
Began listing point-version numbers in this readme.
Added jgi.A_Sample2, an simpler template for a concurrent pipe-filter stage.
Added jgi.MakeChimeras, a tool for making chimeric PacBio reads from input non-chimeric reads.  Also, makechimeras.sh.  Requested by Esther Singer.
Added support for normalized binning to CoveragePileup.  Requested by Vasanth Singan.
v33.59
Fixed pileup's normalized scaling when dealing with 0-coverage scaffolds.
v33.60
Added driver.FilterReadsByName.java and filterbyname.sh.  Allows inclusion or exclusion of reads by name.
Added midpad flag to RandomReads (allows defining inter-scaffold padding).
v33.61
Added ConcurrentReadInputStreamD, prototype for MPI-version of input stream.
Made Read and all classes that might be attached to reads Serializable.
Added DemuxByName and demuxbyname.sh which allows a single file to be split into multiple files based on read names.
v33.62
Added FilterByCoverage and filterbycoverage.sh to filter assemblies based on contig coverage stats (from Pileup).
Added CovStatsLine, an object representation of Pileup's coverage stats.
Added '#' symbol to coverage stats header.
v33.63
Fixed path in filterbycoverage.sh
v33.64
Added custom scripts driver.MergeCoverageOTU and mergeOTUs.sh for Esther.
Added DecontaminateByNormalization, for automating SAG plate decontamination.
Fixed legacy code that set KmerNormalize to use 8 threads in some cases.
Added "fixquality" for capping quality scores at 41. Requested by Bryce Foster.
Added fasta output to kmercountexact. Requested by Alex Copeland.
Added kmer histogram to kmercountexact (2-column and 3-column). Requested by Alex Copeland.
Added multiple memory-related and output formatting flags to kmercountexact.
Made KmerNode a subclass of AbstractKmerTable.
Improved Data's "unloadall" to also clear scaffold-related data.
Removed obsolete class CoverageArray1.
v33.65
Reduced preallocated memory in kmercountexact to avoid a crash on high memory machines.  Also reduced total number of threads.
v33.66
"CountKmersExact.java" renamed to "KmerCountExact.java".
kmercountexact now writes histogram and kmer dump simultaneously in seperate threads.
kmercountexact.sh now specifies both -Xms and -Xmx.
CountKmersExact will no longer run out of memory if -Xms is not specified; instead, it will preallocate a smaller table.
v33.67
Messed with MDA amp in RandomReads a bit.
Added parser "ztd" ("zipthreaddivisor") flag.  Defaults to 2 for removehuman.sh.
Added BBMerge flags "maq" (minaveragequality) and "mee" (mmaxexpectederrors).  Reads violating these will not be attempted to merge.
Added BBMerge "efilter" flag, to allow disabling of the efilter.  Efilter bans merges of reads that have more than the expected number of errors, based on quality scores.
Closed A_Sample2 I/O streams after completion.  Noted by Jon Rood.
Created SynthMDA, a program to make a synthetic MDA'd single cell genome.  This genome would be used as a reference for RandomReads.
Added Reformat "vpair" or (verifypairing) flag, which allows validation of pair names.  Before, it was just interleaved reads.
Pair name validation will now accept identical names, if the "ain" (allowidenticalnames) flag is set.
Updated reformat.sh, repair.sh, bbsplitpairs.sh with new flags.
Removed FastaReadInputStream_old.java.
Added "forcelength" flag to MakeChimeras.
v33.68
Added "ihist" flag to rqcfilter, default "ihist.txt".  Unless this is set to null, BBMerge will run to generate the insert size histogram after filtering completes.
AbstractKmerTable preallocation is now multithreaded.  Unfortunately, this did not result in a speedup.
Added ByteBuilder-related methods to certain Read output formats.
Added ByteStreamWriter.  This is a threaded writer with low overhead, and is substantially faster than TextStreamWriter (perhaps 2x speed).
Fixed a bug in KmerNode (traversing wrong branch during dump).
All AbstractKmerTable subclasses now dump kmers using bsw/ByteBuilder instead of tsw/StringBuilder.
Added ForceTrimLeft/ForceTrimRight flags to Dedupe (requested by Bryce/Seung-Jin).
v33.69
FilterByCoverage (and thus DecontaminatebyNormalization) now produce a log file indicating which contigs were removed.
FilterByCoverage and DecontaminatebyNormalization can now optionally process coverage before and after normalization, and not remove contigs unless the coverage changes by at least some ratio (default 2).  Enable with "mapraw" and optionally "minratio" flag.
Added ihist to file-list.txt.  TODO:  Verify success.
Reads longer than 200bp are now detected as ASCII-33 regardless of their quality values.  This helps with handling PacBio CCS/ROI data.
Added support in FixPairsAndSingles (repair.sh) for reads with names that do not contain whitespace, but still end with "/1" and "/2".
Added qout flag to RandomReads3.
Refactored TextStreamWriter to be more like ByteStreamWriter.
Added gcformat 0 (no base content info printed) to AssemblyStats2 (stats.sh).
v33.70
Updated RQCFilter and BBQC to bring them closer together and improve some of their defaults.  RQCFilter now has more parameters such as k for filtering and trimming.
RQCFilter now correctly produces the insert size histogram.
v33.71
Fixed a bug in Dedupe preventing overlap detection when 'absorb match' and 'absorb containment' were both disabled.  Noted by Shoudan Liang.
Optimized synthetic MDA procedure.
v33.72
Fixed a bug in SynthMDA.java.  Further tweaked parameters.
Added synthmda.sh.
v33.73
Further tweaked SynthMDA defaults to better match some real data sent to me by Shoudan and Alex.
Fixed a bug in BBDuk's mask mode in which all bases in a masked read were assigned quality 0.  Noted by luc (SeqAnswers).
Fixed a small error in KmerCountExact's preallocation calculation.
Added preallocation to BBDuk/BBDuk2.  Not recommended for BBDuk2 because the tables may need unequal sizes.
Added "restrictleft" and "restrictright" flags to BBDuk (not BBDuk2).  These allow only looking for kmer matches in the leftmost or rightmost X bases.  Requested by lankage (SeqAnswers).
v33.74
Added jgi.Shuffle.java to input a read set and output it in random order.  It can also sort by various things (coordinates, sequence, name, and numericID).
Added CallPeaks, which can call peaks from a histogram.  Requested by Kurt LaButti.
Integrated peak calling into BBNorm and KmerCountExact.
BBNorm now has a "histogramcolumns" flag, so it can produce Jellyfish-compatible output.
Added callpeaks.sh.
v33.75
CallPeaks now calls by raw kmer count rather than unique kmer count.  This better detects higher-order peaks.
Finished CrossContaminate.java and added crosscontaminate.sh.
Added "header" and "headerpound" to pileup.sh, to control header presence and whether they start with "#".
Added "prefix" flag to SynthMDA and RandomReads3, to better track origin of reads during cross-contamination trials.
RQCFilter and BBQC now parse 'usejni' flag; rqcfilter.sh and bbqc.sh default to this being enabled.
Added "uselowerdepth" flag to BBNorm (default true).  Allows normalization by depth of higher or lower read.  Set to false by DecontaminateByNormalization.
v33.76
Fixed a bug in synthmda.sh command line.
Fixed build number not being parsed by SynthMDA.
Added some error handling to CrossContaminate, so it shouldn't hang as a result of missing files.
v33.77
SynthMDA now nullifies reference in memory prior to generating reads.
Parser was not correctly setting the number of compression threads when exactly 1 was requested.
Shuffle is now multithreaded, and CrossContaminate defaults to shufflethreads=3.
Shuffle now removes reads as they are printed, reducing memory usage.
Created shellscript templates for generating and assembling full plates of synth MDA data, and ran successfully.
*SamLine was fixed when generating pnext from clipped reads.  Still needs work; pos1 and pos2 need to be recalculated considering clipping.
BBDuk now tracks #contaminant bases as well as #contaminant reads per scaffold for stats.  Additional flag "columns=5" enables this output.
BBDuk stats are now sorted by #bases, not #reads.
BBDuk counting arrays changed from int to long to handle potential overflow.
v33.78
Modified DemuxByName to handle affixes of variable length (though it's less efficient with multiple lengths).
v33.79
Changed the way "pos" and "pnext" are calculated for paired reads to be consistent.  Bug had been noted with soft-clipped reads by Rob Egan.
Changed LOCAL_ALIGN_TIP_LENGTH from 8 to 1.  Previously, soft-clipping would only occur if at least 8 bases would be clipped; not sure why I did that.
Changed the way "tlen" is calculated to compensate for clipping.
v33.80
Changed default decontaminate minratio from2 to 0 (disabling it) because of false negatives.
Changed default decontaminate mincov from 4 to 5 due to a false negative.
Changed default decontaminate kfilter from 63 to 55 to better reflect Spades defaults.
Fixed a bug in filterbycoverage which was outputting contaminant contigs instead of clean contigs.
Added outd (outdirty) flag to FilterByCoverage.
v33.81
Changed decontaminate normalization target from 100 to 50, and minlength from 0 to 500.
Changed decontaminate minc and minp flags from int to float.
v33.82
Changed cross contaminate probability root from 2 to 3 (increasing amount of lower-level contamination).
Fixed a crash bug in sam file generation caused by the change in the way pos was calculated.
v33.83
Added aecc=f, cecc=f, minprob=0.5, depthpercentile=0.8 flags to DecontaminateByNormalization.  Defaults are as listed.
Dropped mindepth to 3 and maxdepth to target; target default changed to 20.
Changed the way mindepth is handled in normalization; now it is based on the depth of the higher read.
v33.84
Added BBNorm prebits flag for setting prefilter cell size (default 2).
Added Decontaminate filterbits and prefilterbits flags, default 32 and 4.  4 was chosen because MDA data has high error kmer counts.
v33.85
Fixed parsing of decontaminate minc and minp (parsed as ints; should have been floats)
Changed default minc to 3.5.
Change default ratio to 1.2.
v33.86
Changed decontaminate default dp to 0.75.
Changed decontaminate default prebits to 2.
Changed decontaminate default minr (min reads) to 20.  Some tiny (~500bp) low-coverage contigs were getting through.
Changed decontaminate mindepth to 2.
Decontaminate results now prints extra columns for read counts and pre-norm coverage.
v33.87
Added "covminscaf" flag to BBMap and Pileup, to supress output of really short contigs.  Default 0.
Changed CrossContaminate coverage distribution from cubic to geometric.
v33.88
Shuffle removing reads caused incredible slowness; it should have set reads to null.  Fixed.
v33.89
Added HashArrayA, HashForestA, KmerNodeA and updated AbstractKmerTable to allow sets of values per kmer.
Refactored all AbstractKmerTable subclasses.
Added scaffold length tracking to BBDuk (for RPKM).
Added RPKM output to BBDuk (enable with "rpkm" flag).
BBDuk now unloads kmers after finishing processing reads.
v33.90
BBDuk counter arrays are now local per-thread, to prevent cache-thrashing.
Added IntList.toString()
Created Seal class, based on BBDuk with values stored in arrays.
Adjusted auto skip settings of BBDuk (increased size threshold for longer skips).
Added BBDuk skip flag (controls minskip and maxskip).
Fixed a bug in DemuxByName/DecontaminateByNormalization/CrossContaminate: attempt to read directories as files.
v33.91
Fixed a bug in BBDuk related to clearing data too early.  Noted by Brian Foster.
v33.92
Added per-reference-file stats counting to BBDuk/Seal, and "refstats" flag.
Added returnList(boolean) to ConcurrentReadStreamInterface.
Removed an extra listen() call from ConcurrentReadInputStreamD.
Documented "addname" flag for stats.sh.
Implemented restrictleft and restrictright for BBDuk2.
Added "nzo" flag for BBDuk/Seal.
Added sdriscoll's reformatted shellscript help for BBDuk and BBMap.  Thanks!
Added more documentation to bbmap.sh (usequality flag).
Added maq (minaveragequality) flag to BBMap, at request of sdriscoll.
Added rename flag to BBDuk/Seal - renames reads based on what sequences they matched.
Added userefnames flag BBDuk/Seal - the names of reference files are used, rather than scaffold IDs.

v33.93
maxindel flag now allows KMG suffix.
Added "speed" flag to BBDuk/Seal.
Added read processing time to BBDuk/Seal output.
BBDuk "fbm" (findbestmatch) mode is now much faster, using variable rather than fixed-length counters.
Fixed BBDuk2 not working when using the "ref" flag rather than "filterref".
Changed AbstractKmerTable subclass names to *1D and *2D.
Made KmerNode a superclass of KmerNode1D and KmerNode2D and eliminated redundant methods.
Eliminated 2D version of HashForest; it now works with 1D and 2D nodes.
Made HashArray a superclass of HashArray1D and HashArray2D.
Created HashArrayHybrid.
Added slow debugging methods to AbstractKmerTable classes, to verify that values were present after being added.
Fixed bug in KmerNode1D; was never changing its value on 'set'.  Probably only affected Seal.  Seal 1D now appears to produce identical output for prealloc and non-prealloc.
Finished debugging KmerNode2D, KmerForest, HashArray2D, HashArrayHybrid, and Seal.
Added "fbm" and "fum" to Seal.
Seal now defaults to 7 ways.
Adjusted Seal's memory preallocation.
Added -Xms flag to BBMergeGapped BBNorm shellscripts.
v33.94
Added -Xms flag to BBDuk and Seal.
Added qskip flag to BBDuk and Seal (for skipping query kmers).
v33.95
Seal now defaults to HashArrayHybrid rather than HashArrayArray2D
v33.96
Fixed a slowdown in Seal and BBDuk caused by sorting list of ID hits.
v33.97
Wrote driver.CorrelateIdentity and matrixtocolumns.sh for identity correlations between 16S and V4.
Wrote jgi.IdentityMatrix and idmatrix.sh for all-to-all alignment.
Added BandedAligner.alignQuadruple() to check all orientations.
BandedAligner now does not clear the full arrays, only the used portion, which can vary depending on read length.
v33.98
No change - build failure.
v33.99
Changed BandedAligner.PenalizeOffCenter().  Indels were getting double-penalized when they led to length mismatches between query and ref.
Added AlignDouble(), but it looks like AlignQuadruple is the only viable method for calculating full identity when the sequences do not start or stop at the same place.
Added test method to ReadStats to ensure the files are safe to write (ReadStats.testFiles()).
Fixed a bug bqhist output giving read 1 and read 2 same values. Noted by Shoudan/Bryce
Fixed a bug in BBDuk initialization when no kmer input supplied. Noted by Bill A.
Fixed a bug in BBDuk/Seal giving a spurious warning.
Detected race condition in ByteFile2 triggered by closing early.  Not very important.
Added jni path flags to BBDuk shellscript command line.
Wrote FindPrimers and msa.sh to locate primer sites.  Uses MultiStateAligner; outputs in sam format.
Wrote CutPrimers and cutprimers.sh to cut regions flanked by mapped primer locations from sequences, e.g. V4.

TODO: Plot correlation of V4 and 16s.
TODO: Add length into edges of Dedupe output. (Ted)
TODO: Benchmark Seal.  Speed seems inconsistent.
TODO: Locking version of Seal.
TODO: HashArray resize - grow fast up to a limit, then resize to exactly the max allowable.
TODO: Alicia BBMap PacBio slowdown (try an older version...)
TODO: BBMerge rename mode with insert sizes.
TODO: Dump info about Seal kmer copy histogram.
TODO: rpkm for pileup / BBMap.
TODO: qskip (query skip for bbduk/seal)
TODO: Dedupe crash bug. (Kurt)
TODO: CallPeaks minwidth should be a subsumption threshold, not creation threshold.
TODO: CallPeaks should not subsume peaks with valleys in between that are very low.
*TODO: Make TextStreamWriter an abstract superclass.
TODO: BBDuk split mode
TODO: Add option for BBMap to convert U to T. (Asaf Levy)
TODO: Add dedupe support for graphing containments and matches.
TODO: Log normalization.
TODO: Prefilterpasses (prepasses)
TODO: Test forcing msa.scoreNoIndels to always run bidirectionally.
TODO: Message for BBNorm indicating pairing (this is nontrivial)
TODO: Cat and Dog filtering
TODO: Average quality for pileup.sh
TODO: Fix ChromArrayMaker which may skip every other scaffold (for now I have reverted to old, correct version). ***Possibly fixed by disabling interleaving; TODO: Test.
TODO: Modify BBDuk to allow trimming of degenerate bases (with the "cd" or "clonedegenerate" flag).
TODO: Consider changing ConcurrentGenericReadInputStream to put read/base statistics into incrementGenerated(), or at least in a function.
TODO: BBSplit produces alignments to the wrong reference in the output for a specific reference. (Shoudan)
TODO: Change the way Ns are handled in cigar strings, both input and output.
TODO: Add #clipped reads/bases to BBMap output.
TODO: Add method for counting number of clipped bases in a read and unclipped length.
TODO: Orientation statistics for BBMap ihist.
TODO: Clarify documentation of 'reads' flag to note that it means reads OR pairs.
TODO: bs flag does not work with BBWrap (Shoudan).
TODO: Make human removal optional in BBQC, and allow aribtrary reference (such as e.coli) for unit testing.
TODO: Fasta input tries to sometimes keep reading from the file when a limited number of reads is specified.  Gives error message but output is fine.
TODO: 'saa' flag sometimes does not work (Shoudan).
TODO: Kmer transition probabilities for binning.
TODO: One coverage file per scaffold; abort if over X scaffolds. (Andrew Tritt)
TODO: Enable JNI by default for BBMap and Dedupe on Genepool.
TODO: Disable cigar string generation when dumping coverage only (?).  This will disable stats, though.
TODO: Pipethread spawned when decompressing from standard in with an external process.
TODO: FileFormat should test interleaving and quality individually on files rather than relying on a static field.
TODO: Refstats (BBSplit) still reports inflated rates for pairs that don't map to the same reference.  This behavior is difficult to change because it is conflated with BBSPlit's output streams.


v32.
Revised all shellscripts to better detect memory in Linux.  This should massively increase reliability and ease of use.
Added append flag.  Allows appending to output files instead of overwriting.
Append flag now should work with BBWrap, with sam files, and with gzipped files.
All statistics are now stored in longs, rather than ints.
Added statistics tracking of # bases as well as # reads.  Updated human-readable output to show 4 columns.
Split bbmerge into gapped (split kmer) and ungapped (overlap only) versions.  bbmerge.sh calls the ungapped version.
Added "qahist" to bbmap - match/sub/ins/del histogram by quality score.
Fixed "pairlen" flag; it was only being used if greater than the default. (Noted by Harald on seqanswers)
Added insert size median and standard deviation to output stats.  The 'ihist=' flag must be set to enable this, otherwise the data won't be tracked. (Requested by Harald on seqanswers)
Fixed bug in which non-ACGTN IUPAC symbols were not being converted to N. (Noted by Leanne on seqanswers)
Changed shellscripts from DOS to Unix EOL encoding.
Added support for "-h" and "--help" in shellscripts (before it was just in java files).
Created Dedupe2 - faster, and supports 1-cluster-per-file output.
Created Dedupe3 - supports more than 2 affix tables.  Uses slightly more memory.
BBMap now generates "sort" shellscripts even if the output is in bam format.
pileup.sh now prints a coverage summary to standard out.
Added 'split' flag to BBMask.
Fixed bug in randomreads allowing paired reads to come from 'nearby' scaffolds.
Documented randomreads.sh.
Added gaussian insert size distribution to randomreads.
Fixed a bug in calcmem.sh that prevented requesting memory that Linux considered 'cached'.
TODO: Penalize score of sites with errors near read tips, and long deletions.
Added "Median_fold" column to pileup.  You need to set 'bitset=
Changed default quality-filtering mode to average probability rather than average quality score.
Default number of threads now takes the environment variable NSLOTS into consideration.  However, because Mendel nodes have hpyerthreading enabled, if NSLOTS>8 and (# processors)==NSLOTS*2, then #processors will be used instead.  So it is still recommended that you set threads manually if you don't have exclusive access to a node.
Fixed bbmerge, which was crashing on fasta input.
Fixed gaussian insert size distribution in randomreads (it was causing a crash).
Enabled unpigz support in Windows (decompression only).
TODO:  BBNorm needs in1/in2/out1/out2 support.
Added mingc and maxgc to reformat.
Added 'passes' flag to BBQC and reduced default passes to 1 if normalization is disabled.
Swapped FileFormat's method signature "allowFileRead" and "allowSubprocess" parms for some functions, as they were inconsistent.  This may have unknown effects.
TODO: unclear if fasta files are currently checked for interleaving.  Method added to "FASTQ".
TODO: FileFormat should perhaps test for quality format and interleaving.
Fixed reversed variables in "machineout" stats for %mapped and %unambiguous.  Found by Michael Barton.
Added "testformat.sh".
Fixed dedupe "csf" output to work even when no other outputs specified.
Fixed dedupe erroneous assumption that "bandwidth" had not been custom-specified.
Changed MakeLengthHistogram (readlength.sh) default behavior to place reads in lower bins rather than closest bins.  Toggle with "round" flag.
Added "repair" flag to SplitPairsAndSingles.  Created "repair.sh".
Fixed a bug in which tabs were not allowed in fasta headers.
Improved BBMerge: default minqo 7->8, made margin a parameter, added 'strict' macro that reduces false positive rate.
Added "samestrand" flag to RandomReads.
Fixed a dedupe bug with "pto" and paired reads; read2 was not getting a UnitID.
Fixed a bug in which the BBMap stats for insertion rate was sometimes higher than the true value.
Fixed bugs in BBMerge; increased speed slightly.
Created grademerge.sh to grade merged reads.
Added 'variance' flag to randomreads; used to make qualities less uniform between reads.
BBDuk now has overwrite=true by default.
calcmem.sh now sets -Xmx and -Xms from each other if only one was specified.
Fixed bug with "ambig=all" and "stoptag" flags being used together.  Found by WhatSoEver (seqanswers).
Added 'findbestmatch'/'fbm' flag to BBDuk; reports the reference sequence sharing the greatest number of kmers with the read.
Shellscripts no longer try to calculate memory before displaying help (noted by Kjiersten Fagnan).
-ea and -da are now valid parameters for all shellscripts.
Improved documentation of Dedupe.
Added "loose" and "vloose" modes to BBMerge.
Added novel-kmer-filtering to BBMerge - bans merged reads that create a novel kmer.  Does not seem to help.
Added entropy-detection to BBMerge - minimum allowed overlap is determined by entropy rather than a constant.  Moderate improvement.
Fixed bug causing "repair.sh" script to not work.  Noted by SES (seqanswers).
Added "fast" mode to BBMerge.
Fixed a rounding problem in RandomReads that caused gaussian distribution to have 2x frequency of intended reads at exactly insert size of double read length.
Added exponential decay insert size distribution to RandomReads, for use in LMP libraries.
TODO: Track different paired read orientation rates (innie, outie, same direction, etc) with BBMap.
Added sssr (secondarysitescoreratio) and ssao (secondarysiteasambiguousonly) flags.  Response to WhatSoEver (seqanswers).
Ambiguously-mapped reads that print a primary site now print a minimum of 1 secondary site, and all sites with the same score as the top secondary site.
Improved error message for paired reads with unequal number of read 1 vs read 2.  Response to Salvatore (seqanswers).
Updated bbcountunique.sh help message.
Changed AddAdapters default to "arc=f" (no reverse-complement adapters).  Added "addpaired" flag (adds adapter to same location of both reads).
Added BBDuk/BBDuk2 "tbo" (trimbyoverlap) flag.  Vastly reduces false-negatives with no increase in false-positives.
Adding "fragadapter" flag to RandomReads.  Also added ability to handle multiple different adapters for both read 1 and read 2.  Adapters are added to paired reads with insert size shorter than read length.
Added "ordered" flag to BBDuk/BBDuk2.
Added "tpe" (trimpairsevenly) flag to BBDuk/BBDuk2.  This only works in conjunction with kimer-trimming to the right.  Slightly decreases false negatives and doubles false positives.
Updated rqcfilter and bbqc with 'tbo' and 'tpe' flags.
TODO: Migrate RQCFilter to BBDuk2.
Improved addadapters to better handle reads annotated by renamereads.
BBMap's fillLimited routine is now affected by 'sssr' flag, if secondary sites are enabled.  This will make things slightly slower when secondary sites are enabled, if sssr uses a low value (default is 0.95).
statswrapper now allows comma-delimited files.
Added standard deviation to BBMerge (requested by Bryce F).
Added "tbo" (trimbyoverlap) flag to BBMerge, as an alternative to joining.
Updated help for 'ambig' in bbmap.sh to remove the obsolete information that 'ambig=all' did not support sam output.
Updated BBMapSkimmer and its shellscript to default to 'ambig=all', which is its intended mode.
BBDuk no longer defaults to "out=stdout.fq" because that was incredibly annoying.  Now it defaults to "out=null".
Changed BBDuk default mink from 4 to 6.
Changed BBDuk, Reformat, SplitPairsAndSingles default trimq from 4 to 6.
Added "ftr"/"ftl" flags to BBDuk.
Added "bbmapskimmer" to the list of options parsed by BBWrap.  (Noted by JJ Chai)
Corrected documentation of idtag and stoptag - both default to false, not true. (Noted by JJ Chai)
Added "mappedonly" flag to reformat. (Requested by Kristen T)
Added "rmn" (requirematchingnames) flag to Dedupe.  Requested by Alex Copeland.
Added ehist, indelhist, idhist, gchist, lhist flags to BBMap, BBDuk, and Reformat.
Added removesmartbell.sh wrapper for pacbio.RemoveAdapters2.
Fixed instance in KmerCoverage where input stream was being started twice.  Noted by Alicia Clum.
Added "ngn" (NumberGraphNodes) flag to dedupe; default true.  Allows toggling of labelling graph nodes with read number or read name.
"slow" flag now disables a heuristic that skipped mapping reads containing only kmers that are highly overrepresented in the reference.  Problem noted by Shoudan Liang.
Added MergeBarcodes and mergebarcodes.sh
Identity is now calculated neutrally by default.
Added "qin" and "qout" documentation to bbnorm shellscripts. Noted by muol (seqanswers).
Changed qhist to ouput additional columns - both linear averages and logrithmic averages.
Added mode to BBMerge output.
Added mode, min, max, median, and standard deviation to ReadLength output.  The mode and std dev are affected by bin size, so will only be exactly correct when bin size is 1.
Added "nzo" (nonzeroonly) flag to ReadLength.
Created "A_Sample", a template for programs that input reads, perform some function, and output reads.
BBNorm now works correctly with dual input and output files.  Noted by Olaf (seqanswers).
Added mode to BBMap insert size statistics.
Added CorrelateBarcodes and filterbarcodes.sh, for analyzing and filtering reads by barcode quality.
Added "aqhist" (average quality histogram) to ReadStats - can be used by BBMap, BBDuk, Reformat.


v31.
TODO:  Change pipethreads to redirects (where possible), and hash pipethreads by process, not by filename.
TODO:  Improve scoring function by using gembal distribution and/or accounting for read length.
TextStreamWriter was improperly testing for output format 'other'.  Noted by Brian Foster.
Fixed bug for read stream 2 in RTextOutputStream3.  Found by Brian Foster.
Fixed bug in MateReadsMT creating an unwanted read stream 2.  Found by Brian Foster.
TrimRead.testOptimal() mode added, and made default when quality trimming is performed; old mode can be used with 'otf=f' flag.
Fixed a couple cases where output file format was set to "ordered" even though the process was singlethreaded; this had caused an out-of-memory crash noted by Bill A.
Changed shellscripts of MapPacBio classes to remove "interleaved=false" term.
Reduced Shared.READ_BUFFER_LENGTH from 500 to 200 and Shared.READ_BUFFER_MAX_DATA from 1m to 500k, to reduce ram usage of buffers.
Noticed small bug in trimming; somehow a read had a 'T' with quality 0, which triggered assertion error.  I disabled the assertion but I'm not sure how it happened.
Fixed bug in which pigz was not used to decompress fasta files.
All program message information now defaults to stderr.
Added "ignorebadquality" (ibq) flag for reads with out-of-range quality.
TODO: mask by information content
Added "mtl"/"mintrimlength" flag (default 60).  Reads will not be trimmed shorter than that.
Made 'tuc' (to uppercase) default to true for bbmap, to prevent assertion errors.  Reads MUST be uppercase to match reference.
Added new tool, BBMask.
Reads and SamLines can now be created with null bases.  
SamLines to Read is now faster, skipping colorspace check.
Added deprecated 'SOH' symbol support to FastaInputStream.  This will be replaced with a '>'.  Needed to process NCBI's NT database.
Added "sampad" or "sp" flag to BBMask, to allow masking beyond bounds of mapped reads.
TODO: %reads with ins, del, splice
TODO: #bases mapped/unmapped, avg read length mapped/unmapped
Dedupe now tracks and prints scaffolds that were duplicates with "outd=".  (request by Andrew Tritt)
Updated all shellscripts to support the -h and --help flags. (suggested by Westerman)
RAM detection is now skipped if user supplies -Xmx flag, preventing a false warning. (noted by Westerman)
Created AddAdapters.java.  Capable of adding adapter sequence to a fastq file, and grading the trimmed file for correctness.
Removed some debug code from FileFormat causing a crash on "stdin" with no extension.  Noted by Matt Nolan.
Added BBWrap and bbwrap.sh.  Wraps BBMap to allow multiple input/output files without reloading the reference.
Added support for breaking long fastq reads into shorter reads (maxlength and minlength flags).  Requested by James Han.
Added Pileup support for residual bins smaller than binsize. Flag "ksb", "keepshortbins".  Requested by Kurt LaButti.
Fixed support for breaking long reads; was failing on the last read in the set.  Noted by James Han.
Improved accuracy slightly by better detecting when padding is needed.
Improved verbose output from MSA.
Created TranslateSixFrames, first step toward amino acid mapping.
Improved RandomReads ability to simulate PacBio error profile.
Fixed crash when using BBSplit in PacBio mode. (Noted by Esther Singer)
May have improved ability to read relatively-pathed files if "." is not in $PATH. (nope, seems not)
Fixed crash when using "usequality=f" flag with fasta input reads. (Noted by Esther Singer)
Corrected behaviour of minlength with regards to trimming; it was not always working correctly.
Added "bhist" (base composition histogram) flag.

v30.
Disabled compression/decompression subprocesses when total system threads allowed is less than 3.
Fixed assertion error in calcCorrectness in which SiteScores are not necessarily sorted if AMBIGUOUS_RANDOM=true.  Noted by Brian Foster.
Fixed bug in toLocalAlignment with respect to considering XY as insertions, not subs.  
TODO: XY should be standardized as substitutions.
Added scarf input support. Requested by Alex Copeland.
TODO: Allow sam input with interleaved flag.
TODO: Make pigz a module dependency or script load.
Fixed bug with nodisk mode dropping the name of the first scaffold of every 500MB chunk after the first.  Noted by Vasanth Singan.
Overhaul of I/O channel creation.  Sequence files are now initialized with a FileFormat object which contains information about the format, permission to overwrite, etc.
Increased limit of number of index threads in Windows in nodisk mode (since disk fragmentation is no longer relevant).
Renamed Read.list to sites; added Read.topSite() and Read.numSites(); replaced many instances of things like "r.sites!=null && !r.sites.isEmpty()"
Refactored to put Read and all read-streaming I/O classes in 'stream' package.
Moved kmer hashing and indexing classes to kmer package.
Moved Variation, subclasses, and related classes to var package.
Moved FastaToChrom and ChromToFasta to dna package.
Moved pacbio error correction classes to pacbio package.
Removed stack, stats, primes, and other packages; prefixed all unused pacakges with z_.
TODO: Sites failing Data.isSingleScaffold() test should be clipped, not discarded.
RandomReads3 no longer adds /1 and /2 to paired fastq read names by default (can be enabled with 'addpairnum' flag).
Added "inserttag" flag; adds the insert size to sam output.
Fixed insert size histogram anomaly.  There was a blip at insert==(read1.length+read2.length) because the algorithm used to calculate insert size was different for reads that overlap and reads that don't overlap.
Skimmer now defaults to cigar=true.
Added maxindel1 and maxindel2 (or maxindelsum) flags.
Removed OUTER_DIST_MULT2 because it caused assertion errors when different from OUTER_DIST_MULT; changed OUTER_DIST_MULT from 15 to 14.
Added shellscript for skimmer, bbmapskimmer.sh
TODO:  Document above changes to parameters.



v29.
New version since major refactoring.
Added FRACTION_GENOME_TO_EXCLUDE flag (fgte).  Setting this lower increases sensitivity at expense of speed.  Range is 0-1 and default is around 0.03.
Added setFractionGenometoExclude() to Skimmer index.
LMP librares were not being paired correctly.  Now "rcs=f" may be used to ignore orientation when pairing.  Noted by Kurt LaButti.
Allocating memory to alignment score matrices caused uncaught out-of-memory error on low-memory machines, resulting in a hang.  This is now caught and results in an exit.  Noted by Alicia Clum.
GPINT machines are now detected and restricted to 4 threads max.  This helps prevent out-of-memory errors with PacBio mode.
Fixed sam output bug in which an unmapped read would get pnext of 0 rather than 1 when its mate mapped off the beginning of a scaffold.  Noted by Rob Egan.
Added memory test prior to allocating mapping threads.  Thread count will be reduced if there is not enough memory.  This is to address the issue noted by James Han, in which the PacBio versions would crash after running out of memory on low-memory nodes.
TODO:  Detect and prevent low-memory crashes while loading the index by aborting.
Fixed assertion error caused by strictmaxindel mode (noted by James Han).
Added flag "trd" (trimreaddescriptions) which truncates read names at the first whitespace.
Added "usequality/uq" flag to turn on/off usage of quality information when mapping.  Requested by Rob Egan.
Added "keepbadkeys/kbk" flag to prevent discarding of keys due to low quality.  Requested by Rob Egan.
Fixed crash with very long reads and very small kmers due to exceeding length of various kmer array buffers.
Avg Initial Sites and etc no longer printed for read 2 data.
TODO:  Support for selecting long-mate-pair orientation has been requested by Alex C.
Fixed possible bug in read trimming when the entire read was below the quality threshold.
Fixed trim mode bug: "trim=both" was only trimming the right side.  "qtrim" is also now an alias for "trim".
Fixed bug in ConcurrentGenericReadInputStream causing an incorrect assertion error for input in paired files and read sampling.  Found by Alex Copeland.
Added insert size histogram: ihist=<file>
Added "machineout" flag for machine-readable output stats.
TODO: reads_B1_100000x150bp_0S_0I_0D_0U_0N_interleaved.fq.gz (ecoli) has 0% rescued for read1 and 0.7% rescued for read 2.  After swapping r1 and r2, .664% of r2 is rescued and .001% of r1 is rescued.  Why are they not symmetric?
Added 'slow' flag to bbmap for increased accuracy.  Still in progress.
Added MultiStateAligner11ts to MSA minIdToMinRatio().
Changed the way files are tested for permission to write (moved to Tools).
Fixed various places in which version string was parsed as an integer.
Added test for "help" and "version" flags.
Fixed bug in testing for file existence; noted by Bryce Foster.
Fixed issue with scaffold names not being trimmed on whitespace boundaries when 'trd=t'.  Noted by Rob Egan.
Added pigz (parallel gzip) support, at suggestion of Rob Egan.
Improved support for subprocesses and pipethreads; they are now automatically killed when not needed, even if the I/O stream is not finished.  This allows gunzip/unpigz when a file is being partially read.
Added shellscript test for the hostname 'gpint'; in that case, memory will be capped at 4G per process.
Changed the way cris/ros are shut down.  All must now go through ReadWrite.closeStreams()
TODO:  Force rtis and tsw to go through that too.
TODO:  Add "Job.fname" field.
Made output threads kill processes also.
Modified TrimRead to require minlength parameter.
Fixed a bug with gathering statistics in BBMapPacBioSkimmer (found by Matt Scholz).
Fixed a bug in which reads with match string containing X/Y were not eligible to be semiperfect (Found by Brian Foster).
Fixed a bug related to improving the prior fix; I had inverted an == operator (Found by Brian Foster).
Added SiteScore.fixXY(), a fast method to fix reads that go out-of-bounds during alignment.  Unfinished; score needs to be altered as a result.
Added "pairsonly" or "po" flag.  Enabling it will treat unpaired reads as unmapped, so they will be sent to 'outu' instead of 'outm'.  Suggested by James Han and Alex Copeland.
Added shellscript support for java -Xmx flag (Suggested by James Han).
Changed behavior: with 'quickmatch' enabled secondary sites will now get cigar strings (mostly, not all of them).
"fast" flag now enables quickmatch (50% speedup in e.coli with low-identity reads).  Very minor effect on accuracy.
Fixed bug with overflowing gref due GREFLIMIT2_CUSHION padding.  Found by Alicia Clum.
Fixed bug in which writing the index would use pigz rather than native gzip, allowing reads from scaffolds.txt.gz before the (buffered) writing finished.  Rare race condition.  Found by Brian Foster.
Fixed stdout.fa.gz writing uncompressed via ReadStreamWriter.
Added "allowSubprocess" flag to all constructors of TextFile and TextStreamWriter, and made TextFile 'tryAllExtensions' flag the last param.
allowSubprocess currently defaults to true for ByteFiles and ReadInput/Output Streams.
TODO: TextFile and TextStreamWriter (and maybe others?) may ignore ReadWrite.killProcess().
TODO: RTextOutputStream3 - make allowSubprocess a parameter
TODO: Assert that first symbol of reference fasta is '>' to help detect corrupt fastas.
Improved TextStreamWriter, TextFile, and all ReadStream classes usage of ReadWrite's InputStream/OutputStream creation/destruction methods.
All InputStream and OutputStream creation/destruction now has an allowSubprocesses flag.
Added verbose output to all ReadWrite methods.
Fixed bug in which realigned SiteScores were not given a new perfect/semiperfect status.  Noted by Brian Foster and Will Andreopoulos.


v28.
New version because the new I/O system seems to be stable now.
Re-enabled bam input/output (via samtools subprocess).  Lowered shellscript memory from 85% to 84% to provide space for samtools.
Added "-l" to "#!/bin/bash" at top.  This may make it less likely for the environment to be messed up.  Thanks to Alex Boyd for the tip.
Addressed potential bug in start/stop index padding calculation for scaffolds that began or ended with non-ACGT bases.
Made superclass for Index.
Made superclass for BBMap.
Removed around 5000 lines of code as a result of dereplication into superclasses.
Added MultiStateAligner11ts, which uses arrays for affine transform instead of if blocks.  Changing insertions gave a ~5% speedup; subs gave an immeasurably small speedup.
Found bug in calculation of insert penalties during mapping.  Fixing this bug increases speed but decreases accuracy, so it was modified toward a compromise.


v27.
Added command line to sam file header.
Added "msa=" flag.  You can specify which msa to use by entering the classname.
Added initial banded mode.  Specify "bandwidth=X" or "bandwidthratio=X" accelerate alignment.
Cleaned up argument parsing a bit.
Improved nodisk mode; now does not use the disk at all for indexing.  BBSplitter still uses the disk.
Added "fast" flag, which changes some paramters to make mapping go faster, with slightly lower sensitivity.
Improved error handling; corrupt input files should be more likely to crash with an error message and less likely to hang.  Noted by Alex Copeland.
Improved SAM input, particularly coordinates and cigar-string parsing; this should now be correct but requires an indexed reference.  Of course this information is irrelevant for mapping so this parsing is turned off by default for bbmap.
Increased maximum read speed with ByteFile2, by using 2 threads per file.  May be useful in input-speed limited scenarios, as when reading compressed input on a node with many cores.  Also accelerates sam input.
TODO:  Consider moving THREADS to Shared.
Updated match/cigar flag syntax.
Updated shellscript documentation.
Changed ByteFile2 from array lists to arrays; should reduce overhead.
TODO:  Increase speed of sam input.
TODO:  Increase speed of output, for all formats.
TODO:  Finish ReadStreamWriter.addStringList(), which allows formatting to be done in the host.
In progress:  Moving all MapThread fields to abstract class.
MapThread now passes reverse-complemented bases to functions to prevent replication of this array.
Fixed very rare bug when a non-semiperfect site becomes semiperfect after realignment, but subsequently is no longer highest-ranked.
strictmaxindel can now be assigned a number (e.g. stricmaxindel=5).
If a fasta read is broken into pieces, now all pieces will recieve the _# suffix in their name.  Previously, the first piece was exempt.
TODO:  Consider changing SamLine.rname to a String and seq, qual to byte[].
Changed SamLine.seq, qual to byte[].  Now stored in original read order and only reversed for minus strand during I/O.
Added sortscaffolds flag (requested by Vasanth Singan).
Fixed XS tag bug; in some cases read 2 was getting opposite flag (noted by Vasanth Singan).
Fixed bug when reading sam files without qualities (noted by Brian Foster).
Fixed bug where absent cigar strings were printed as "null" instead of "*" as a result of recent changes to sam I/O (noted by Vasanth Singan).
Found error when a read goes off the beginning of a block.  Ref padding seems to be absent, because Ns were replaced by random sequence.  Cause is unknown; cannot replicate.
Fixed Block.getHitList(int, int).
Changed calcAffineScore() to require base array for information when throwing exceptions.
Changed generated bamscript to unload samtools module before loading samtools/0.1.19.
sam file idflag and stopflag are both now faster, particularly for perfect mappings.  But both default to off because they are still slow nonetheless.
Fixed bug in BBIndex in which a site was considered perfect because all bases matched the reference, but some of the bases were N.  Canonically, reads with Ns can never be perfect even if the ref has Ns in the same locations.
Fixed above bug again because it was not fully fixed:  CHECKSITES was allowing a read to be classified as perfect even if it contained an N.
Increased sam read speed by ~2x; 30MB/s to 66MB/s
Increased sam write speed from ~18MB/s to ~32MB/s on my 4-core computer (during mapping), with mapping at peak 42MB/s with out=null.  Standalone (no mapping) sam output seems to run at 51MB/s but it's hard to tell.
Increased fasta write from 118MB/s to 140 MB/s
Increased fastq write from 70MB/s to 100MB/s
Increased fastq read from 120MB/s (I think) to 296MB/s (663 megabytes/sec!) with 2 threads or 166MB/s with 1 thread
Some of these speed increases come from writing byte[] into char[] buffer held in a ThreadLocal, instead of turning them into Strings or appending them byte-by-byte.
All of these speed optimizations caused a few I/O bugs that temporarily affected some users between Oct 1 and Oct 4, 2013.  Sorry!
Flipped XS tag from + to - or vice versa.  I seem to have misinterpreted the Cufflinks documentation (noted by Vasanth Singan).
Fixed bug in which (as a result of speed optimizations) reads outside scaffold boundaries, in sam 1.3 format, were not getting clipped (Noted by Brian Foster).
Changed default behavior of all shellscripts to run with -Xmx4g if maximum memory cannot be detected (typically, because ulimit=infinity).  Was 31.  Unfortunately things will break either way.
Fixed off-by-1 error in sam TLEN calculation; also simplified it to give sign based on leftmost POS and always give a plus and minus even when POS is equal.
Added sam NH tag (when ambig=all).
Disabled sam XM tag because the bowtie documentation and output do not make any sense.
Changed sam MD and NM tags to account for 'N' symbol in cigar strings.
Made sam SM tag score compatible with mapping score.
Fixed bug in SamLine when cigar=f (null pointer when parsing match string). (Found by Vasanth Singan)
Fixed bug in BBMapThread* when local=true and ambiguous=toss (null pointer to read.list). (Found by Alexander Spunde)
Changed synthetic read naming and parsing (parsecustom flag) to use " /1" and " /2" at the end of paired read names. (Requested by Kurt LaButti)
Increased fastq write to 200MB/s (590 megabytes/s)
Increased fasta write to 212MB/s (624 megabytes/s measured by fastq input)
Increased sam write to 167MB/s (492 megabytes/s measured by fastq input)
Increased bread write to 196MB/s (579 megabytes/s measured by fastq input)
bf2 (multithreaded input) is now enabled by default on systems with >4 cores, or in ReformatReads always.
Fixed RTextOutputStream3.finishedSuccessfully() returning false when output was in 2 files.
Changed output streams to unbuffered.  No notable speed increase.
Fixed bug in ByteFile2 in which reads would be recycled when end of file was hit (found by Brian Foster, Bryce Foster, and Kecia Duffy).


v26.
Fixed crash from consecutive newlines in ByteFile.
Made SiteScore clonable/copyable.
Removed @RG line from headers.  It implies that reads should be annotated with addition fields based on the RG line information.
Changed sam flags (at advice of Joel Martin).  Now single-ended reads will never have flags 0x2, 0x40, or 0x80 set.
Added correct insert size average to output stats, in place of old inner distance and mapping length.
Fixed crash when detecting length of SamLines with no cigar string.  (Found by Shayna Stein)
Added flag "keepnames" which keeps the read names unchanged when writing in sam format.  Normally, a trailing "/1", "/2", " 1", or " 2" are stripped off, and if read 2's name differs from read 1's name, read 1's name is used for both.  This is to remain spec-compliant with the sam format.  However, in some cases (such as grading synthetic reads tagged with the correct mapping location) it is useful to retain the original name of each read.
Added local alignment option, "local".  Translates global alignments into a local alignments using the same affine transform (and soft-clips ends).
Changed killbadpairs default to false.  Now by default improperly paired reads are allowed.
Merged TranslateColorspaceRead versions into a single class.
Added interleaved input and output for bread format.  May be useful for error correction pipeline.
TODO:  Mode where reads are mapped to multiple scaffolds, but are mapped at most one time per scaffold.  I.e., remove all but top site per scaffold (and forbid self-mapping).
Fixed yet another instance of negative coordinates appearing in an unmapped read, which the new version of samtools can't handle.
Fixed bug in counting ambiguous reads; was improperly including in statistics reads that were ambiguous but had a score lower than minratio.
Fixed rare crash found related to realignment of reads with ambiguous mappings (found by Rob Egan).
Unified many of the differences between the MapThread variants, and added a new self-checking function (checkTopSite) to ensure a Read is self-consistent.
Added some bitflag fetch functions to SamLine and fixed 'pairedOnSameChrom()' which was not handling the '=' symbol.
TODO: Make GENERATE_BASE_SCORES_FROM_QUALITY a parameter, default false in BBMapPacBio and true elsewhere. (I verified this should work fine)
TODO: Make GENERATE_KEY_SCORES_FROM_QUALITY a parameter, default true (probably even in BBMapPacBio). (I verified this should work fine)
Updated LongM (merged with LongM from Dedupe).
Fixed bug in SamLine in which clipped leading indels were not considered, causing potential negative coordinates.  (Found by Brian Foster)
TODO: Match strings like NNNNNNDDDDDNNNNNmmmmmmmmmmmmmmmmm...mmmmmmm should never exist in the first place.  Why did that happen?
Added "strictmaxindel" flag (default: strictmaxindel=f).  Attempts to kill mappings in which there is a single indel event longer than the "maxindel" setting.  Requested by James Han.
TODO: Ensure strictmaxindel works in all situations, including rescued paired ends and recursively regenerated padded match strings.
TODO: Redo msa to be strictly subtractive.  Start with score=100*bases, then use e.g. 0 for match, -1 for del, -370 for sub, -100 for N, etc.  No need for negative values.
Changed TIMEBITS in MultiStateAligner9PacBio from 10 to 9 to address a score underflow assertion error found by Alicia Clum.  The underflow occuerd around length 5240; new limit should be around 10480.
TODO: Alicia found an error of exceeding gref bounds.
Fixed race condition in TextStreamWriter.
Improved functionality of splitter.  Now you can index once and map subsequently using "basename" without specifying "ref=" every single time.
"Reads Used" in output now dispays the number of reads used.  Before, for paired reads, it would display the number of pairs (half as many).
Added bases used to reads used at Kurt's request.
Improved bam script generation.  Now correctly sets samtools memory based on detected memory, and warns user that crashes may be memory-related.
Fixed an obsolete assertion in SamLine found by Alicia.
Added XS tag option ("xstag=t") for Cufflinks; the need for this was noted by requested by Vasanth Singan.
Added 'N' cigar operation for deletions longer than X bases (intronlen=X).  Also needed by Cufflinks.
Secondary alignments now get "*" for bases and qualities, as recommended by the SAM spec.  This saves space, but may cause problems when converting sam into other formats.
Fixed bug that caused interleaved=true to override in2.  Now if you set in and in2, interleaved input will be disabled.  (noted by Andrew Tritt).
Fixed some low-level bugs in I/O streams.  When shutting down streams I was waiting until !Thread.isAlive() rather than Thread.getState()==Thread.State.TERMINATED, which caused a race condition (since a thread is not alive before it starts execution).
Added debugging file with random name written to /ref/ directory.  This should help debugging if somewhere deep in a pipeline multiple processes try to index at the same location simultaneously.  Suggested by Bryce Foster.
Fixed log file generation causing a crash if the /ref/ directory did not exist, found by Vasanth Singan.  Also logging is now disabled by default but enabled if you set "log=t".
Input sequence data will now translate '.' and '-' to 'N' automatically, as some fasta databases appear to use '.' instead of 'N'.  (Thanks to Kecia Duffy and James Han)
Added capability to convert lowercase reads to upper case (crash on lowercase noted by Vasanth Singan).


v25.
Increased BBMapPacBio max read length to 6000, and BBMapPacBioSkimmer to 4000.
Fixed bugs in padding calculations during match string generation.
Improved some assertion error output.
Added flag "maxsites" for max alignments to print.
Added match field to sitescore.
Made untrim() affect sitescores as well.
Decreased read array buffer from 500 to 20 in MapPacBio.
TODO:  stitcher for super long reads.
TODO:  wrapper for split reference mapping and merging.
Improved fillAndScoreLimited to return additional information.
Added flag "secondary" to print secondary alignments.  Does not yet ensure that all secondary alignments will get cigar strings, but most do.
Added flag "quickmatch" to generate match strings for SiteScores during slow align.  Speeds up the overall process somewhat (at least on my PC; have not tested it on cluster).
Improved pruning during slow align by dynamically increasing msa limit.
Addressed a bug in which reads sometimes have additional sites aligned to the same coordinates as the primary site.  The bug can still occur (typically during match generation or as a result of padding), but is detected and corrected during runtime.
Tracked down and fixed a bug relating to negative coordinates in sam output for unmapped reads paired with reads mapped off the beginning of a scaffold, with help from Rob Egan.
Disabled frowny-face warning message which had caused some confusion.
TODO:  Add verification of match strings on site scores.
Made superclass for MSA.  This will allow merging of redundant code over the various BBMap versions.
Fixed a crash-hang out-of-memory error caused by initialization order.  Now crashes cleanly and terminates.  Found by James Han.
Fixed bug in output related to detecting cigar string length under sam 1.4 specification (found by Rob Egan).
Added flag "killbadpairs"/"kbp".
Added flag "fakequality" for fasta.
Permanently fixed bugs related to unexpected short match strings caused by error messages.
Increased speed of dynamic program phase when dealing with lots of Ns.
TODO:  In-line generation of short match string when printing a read, rather than mutating the read. (mutation is now temporary)
Added flag, "stoptag".  Allows generation of SAM tag YS:i:<read stop location>
Added flag, "idtag".  Allows generation of SAM tag YI:f:<percent identity>

v24.
Fixed bug that slightly reduced accuracy for reads with exactly 1 mismatch.  They were always skipping slow align, sometimes preventing ambiguous reads from being detected.
Increased speed of MakeRocCurve (for automatic grading of sam files from synthetic reads).  Had used 1 pass per quality level; now it uses only 1 pass total.
Increased accuracy of processing reads and contigs with ambiguous bases (in mapping phase).
Adjusted clearzones to use gradient functions and asymptotes rather than step functions.  Reduces false positives and increases true positives, especially near the old step cutoffs.
Fixed trimSitesBelowCutoff assertion that failed for paired reads.
Added single scaffold toggle to RandomReads.  Default 'singlescaffold=true'; forces reads to come from a single scaffold).  This can cause non-termination if no scaffolds are long enough, and may bias against shorter scaffolds.
Added min scaffold overlap to RandomReads.  Default 'overlap=1'; forces reads to overlap a scaffold at least this much.  This can cause non-termination if no scaffolds are long enough, and may bias against shorter scaffolds.
Fixed setPerfect().  Previously, reads with 'N' overlapping 'N' in the reference could be considered perfect matches, but no reads containing 'N' should ever be considered a perfect mapping to anything.
Formalized definition of semiperfect to require read having no ambiguous bases, and fixed "isSemiperfect()" function accordingly.
Shortened and clarified executable names.
Fixed soft-clipped read start position calculation (mainly relevant to grading).
Prevented reads from being double-counted when grading, when a program gives multiple primary alignments for a read.
Fixed a bug in splitter initialization.
Added "ambiguous2".  Reads that map to multiple references can now be written to distinct files (prefixed by "AMBIGUOUS_") or thrown away, independantly of whether they are ambiguous in the normal sense (which includes ambiguous within a single reference).
Added statistics tracking per reference and per scaffold.  Enable with "scafstats=<file>" or "refstats=<file>".
"ambiguous" may now be shortened to "ambig" on the command line.
"true" and "false" may now be shortened to t, 1, or f, 0.  If omitted entirely, "true" is assumed; e.g. "overwrite" is equivalent to "overwrite=true".
Added stderr as a vaild output destination specified from the command line.
BBSplitter now has a flag, "mapmode"; can be set to normal, accurate, pacbio, or pacbioskimmer.
Fixed issue where stuff was being written to stdout instead of stderr and ended up in SAM files (found by Brian Foster).
TODO: Add secondary alignments.
TODO: Unlimited length reads.
TODO: Protein mapping.
TODO: Soft clipping in both bbmap and GradeSamFile.  Should universally adjust coords by soft-clip amount when reported in SAM format.
Fixed assertion error concerning reads containing Ns marked as perfect, when aligned to reference Ns (found by Rob Egan).
Fixed potential null-pointer error in "showprogress" flag.

v23.
Created BBSplitter wrapper for BBMap that allows merging any number references together and splitting the output into different streams.
Added support for ambiguous=random with paired reads (before it was limited to unpaired).
TODO: Iterative anchored alignment for very long reads, with a full master gref.
TODO: untrim=c/m/s/n/r
TODO: mode=vfast/veryfast: k=14 minratio=0.8 minhits=2 maxindel=20
TODO: mode=fast: k=13 minratio=0.7 minhits=2 maxindel=200
TODO: mode=normal: k=13 minratio=0.56 minhits=1 maxindel=16000
TODO: mode=slow/accurate: BBMapi
TODO: mode=pacbio: BBMapPacBio k=12
TODO: mode=perfect
TODO: mode=semiperfect
TODO: mode=rnaseq
TODO: Put untrim in caclStatistics section
TODO: Test with MEGAN.
Finished new random read generator.  Much faster, and solves coordinate problem with multiple indels.
Improved error message on read parsing failures.
TODO: Insert size histogram
TODO: "outp=", output for reads that mapped paired
TODO: "outs=", output for reads that mapped singly
Corrected assertion in "isSingleScaffold()"
Fixed a rare bug preventing recursive realignment when ambiguous=random (found by Brian Foster)
Added samversion/samv flag.  Set to 1.3 for cigar strings with 'M' or 1.4 for cigar strings with '=' and 'X'.  Default is 1.3.
Added enforcement of thread limit when indexing.
Added internal autodetection of gpint machines.  Set default threadcount for gpints at 2.
Improved ability to map with maxindel=0
Added XM:i:<N> optional SAM flag because some programs seem to demand it.  Like all extra flags, this is omitted if the read is not mapped.  Otherwise, it is set to 1 for unambiguously mapped reads, and 2 or more for ambiguously mapped reads.  The number can range as high as the total number of equal-scoring sites, but this is not guaranteed unless the "ambiguous=random" flag is used.
Fixed bug in autodetection of paired ends, found by Rob Egan.



v22.
Added match histogram support.
Added quality histogram support.
Added interleaving support to random read generator.
Added ability to disable pair rescue ("rescue=false" flag), which can speed things up in some cases.
Disabled dynamic-programming slow alignment phase when no indels are allowed.
Accelerated rescue in perfect and semiperfect mode.
Vastly accelerated paired mapping against references with a very low expected mapping rate.
Fixed crash in rescue caused by reads without quality strings (e.g. paired fasta files). (found by Brian Foster)


v21.
If reference specified is same as already-processed reference, the old index will not be deleted.
Added BBMap memory usage estimator to assembly statistics tool:  java -Xmx120m jgi.AssemblyStats2 <fasta file> k=<kmer size for BBMap>
Added support for multiple output read streams: all reads (set by out=), mapped reads (set by outm=), and unmapped reads (set by outu=).  They can be in different formats and any combination can be used at once.  You can set pair output to secondary files with out2, outm2, and outu2.
Changed definition of "out=".  You can no longer specify split output streams implicitly by using a "#" in the filename; it must be explicit.  the "#" wildcard is still allowed for input streams.
Fixed a bug with sam input not working. (found by Brian Foster)
Added additional interleaved autodetection pattern for reads named "xxxxx 1:xxxx" and "xxxxx 2:xxxx"
Fixed a bug with soft-clipped deletions causing an incorrect cigar length. (found by Brian Foster)
Fixed a bug with parsing of negative numbers in byte arrays.
TODO: Found a new situation in which poly-N reads preferentially map to poly-N reference (probably tip search?)
Fixed a bug in which paired reads occasionally are incorrectly considered non-semiperfect. (found by Brian Foster)
Added more assertion tests for perfection/imperfection status.
Added blacklist support.  This allows selection of output stream based on the name of the scaffold to which a read maps.
Created Blacklist class, allowing creation of blacklists and whitelists.
Added outb (aka outblacklist) and outb2 streams, to output reads that mapped to blacklisted scaffolds. 
Added flag "outputblacklisted=<true/false>" which contols whether blacklisted reads are printed to the "out=" stream.  Default is true.
Added support for streaming references. e.g. "cat ref1.fa ref2.fa | java BBMap ref=stdin.fa"
Updated and reorganized this readme.
Removed a dependency on Java 7 libraries (so that the code runs in Java 6).
Added per-read error rate histogram.  Enable with qhist=<filename>
TODO: generate standard deviation.
Added per-base-position M/S/D/I/N rate tracking.  Enable with mhist=<filename>
Added quality trimming.  Reads may be trimmed prior to mapping, and optionally untrimmed after mapping, so that no data is lost.  Trimmed bases are reported as soft-clipped in this case.
Trimming will extend until at least 2 consecutive bases have a quality greater than trimq (default 5).
Added flags:  trim=<left/right/both/false>, trimq=<5>, untrim=<true/false>
TODO:  Correct insert size in realtime for trim length.
TODO:  Consider adding a TrimRead pointer to reads, rather than using obj.
TODO:  Consider extending match string as 'M' rather than 'C' as long as clipped bases match.
Found and made safe some instances where reads could be trimmed to less than kmer length.
Found and fixed instance where rescue was attempted for length-zero reads.
Fixed an instance where perfect reads were not marked perfect (while making match string).


v20.1 (not differentiated from v20 since the differences are minor)
Fixed a minor, longstanding bug that prevented minus-strand alignment of rads that only had a single valid key (due to low complexity or low quality).
Increased accuracy of perfectmode and semiperfectmode, by allowing mapping of reads with only one valid key, without loss of speed.  They still don't quite match normal mode since they use fewer keys.
Added detection of and error messages for reads that are too long to map. 
Improved shell script usage information.


v20.
Made all MapThreads subclasses of MapThread, eliminating duplicate code.
Any exception thrown by a MapThread will now be detected, allowing the process to complete normally without hanging.
Exceptions (e.g. OutOfMemory) when loading reference genome are now detected, typically causing a crash exit instead of a hang.
Exceptions (e.g. OutOfMemory) when generating index are now detected, causing a crash exit instead of a hang.
Exceptions in output stream (RTextOutputStream) subthreads are now detected, throwing an exception.
Added support for soft clipping.  All reads that go off the ends of scaffolds will be soft-clipped when output to SAM format. (The necessity of this was noted by Rob Egan, as negative scaffold indices can cause software such as samtools to crash)


v19.
Added support for leading FASTA comments (denoted by semicolon).
Fixed potential problem in FASTA read input stream with very long reads.
Recognizes additional FASTA file extensions: .seq, .fna, .ffn, .frn, .fsa, .fas
Disabled gzip subprocesses to circumvent a bug in UGE: Forking can cause a program to be terminated.  Gzip is still supported.
Slightly reduced memory allocation in shellscript.
Ported "Analyze Index" improvement over to all versions (except v5).
Added flags: fastaminread, showprogress
Fixed problem noted by Rob Egan in which paired-end reads containing mostly 'N' could be rescued by aligning to the poly-N section off the end of a contig.
Fixed: Synthetic read headers were being improperly parsed by new FASTQ input stream.
Made a new, faster, more correct version of "isSemiperfect".
Added "semiperfect" test for reads changed during findDeletions.
Identified locations in "scoreNoIndels" where call 'N' == ref 'N' is considered a match.  Does not seem to cause problems.
Noted that SAM flag 0x40 and 0x80 definitions differ from my usage.


v18.
Fastq read input speed doubled.
Fasta read input speed increased 50%.
Increased speed of "Analyze Index" by a factor of 3+ (just for BBMap so far; have not yet ported change over to other versions).
Fixed an array out-of-bounds bug found by Alicia Clum.
Added bam output option (relies on Samtools being installed).
Allows gzip subprocesses, which can sometimes improve gzipping and gunzipping speed over Java's implementation (will be used automatically if gzip is installed).  This can be disabled with with the flags "usegzip=false" and "usegunzip=false".
Started a 32-bit mode which allows 4GB per block instead of 2GB, for a slight memory savings (not finished yet).
Added nondeterministic random read sampling option.
Added flags: minscaf, startpad, stoppad, samplerate, sampleseed, kfilter, usegzip, usegunzip


v17.
Changed the way error rate statistics are displayed.  All now use match string length as denominator.
Identified error in random read generator regarding multiple insertions.  It will be hard to fix but does not matter much.
Found out-of-bounds error when filling gref.  Fixed (but maybe not everywhere...).
Added random mapping for ambiguous reads.
Changed index from 2d array to single array (saves a lot of memory).
Increased speed by ~10%.
Improved index generation and loading speed (typically more than doubled).
Changed chrom format to gzipped.
Added "nodisk" flag; index is not written to disk.
Fixed a rare out-of-bounds error.
Increased speed of perfect read mapping.
Fixed rare human PAR bug.


v16.  Changes since last version:
Supports unlimited number of unscaffolded contigs.
Supports piping in and out.  Set "out=stdout.sam" and "in=stdin.fq" to pipe in a fastq file and pipe out a sam file (other extensions are also supported).
Ambiguously named files (without proper extensions) will be autodetected as fasta or fastq (though I suggest not relying on that).
Added additional flags (described in parameters section): minapproxhits, padding, tipsearch, maxindel.
minapproxhits has a huge impact on speed.  Going from 1 to 2 will typically at least double the speed (on a large genome) at some cost to accuracy.


v15.  Changes since last version:
Contig names are retained for output.
SAM header @SQ tags fixed.
SAM header @PG tag added.
An out-of-bounds error was fixed.
An error related to short match strings was found and possibly handled.
All versions now give full statistics related to %matches, %substitutions, %deletions, and %insertions (unless match string generation is disabled).
Increased speed and accuracy for tiny (<20MB) genomes.
Added dynamic detection of scaffold sizes to better partition index, reducing memory in some cases.
Added command-line specification of kmer length.
Added more command line flags and described them in this readme.
Allowed overwriting of existing indices, for ease of use (only when overwrite=true).  For efficiency you should still only specify "ref=" the first time you map to a particular reference, and just specify the build number subsequently.
