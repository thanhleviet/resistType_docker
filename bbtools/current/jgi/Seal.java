package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;

import kmer.AbstractKmerTable;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

import align2.IntList;
import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import dna.AminoAcid;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;

/**
 * SEAL: Sequence Expression AnaLyzer
 * Derived from BBDuk.
 * Allows multiple values stored per kmer.
 * Intended for RNA-seq or other reads-per-sequence quantification.
 * @author Brian Bushnell
 * @date November 10, 2014
 *
 */
public class Seal {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		//Create a new Seal instance
		Seal pup=new Seal(args);
		
		///And run it
		pup.process();
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx20g -cp <path> jgi.Seal in=<input file> out=<output file> ref=<contaminant files>");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("in2=<file>         \tUse this if 2nd read of pairs are in a different file.");
		outstream.println("out=<file>         \t(outmatch) 'out=stdout' will pipe to standard out.");
		outstream.println("out2=<file>        \t(outmatch2) Use this to write 2nd read of pairs to a different file.");
		outstream.println("outu=<file>        \t(outunmatched) Write unmatched reads here.");
		outstream.println("outu2=<file>       \t(outunmatched2) Use this to write 2nd read of pairs to a different file.");
		outstream.println("stats=<file>       \tWrite statistics about which contaminants were detected.");
		outstream.println("rpkm=<file>        \tWrite coverage and RPKM/FPKM info.");
		outstream.println("");
		outstream.println("threads=auto       \t(t) Set number of threads to use; default is number of logical processors.");
		outstream.println("overwrite=t        \t(ow) Set to false to force the program to abort rather than overwrite an existing file.");
		outstream.println("showspeed=t        \t(ss) Set to 'f' to suppress display of processing speed.");
		outstream.println("interleaved=auto   \t(int) If true, forces fastq input to be paired and interleaved.");
		outstream.println("k=31               \tKmer length used for finding contaminants.  Contaminants shorter than k will not be found.");
		outstream.println("maskmiddle=t       \t(mm) Treat the middle base of a kmer as a wildcard.");
		outstream.println("minkmerhits=0      \t(mh) Reads with more than this many contaminant kmers will be discarded.");
		outstream.println("minavgquality=0    \t(maq) Reads with average quality (before trimming) below this will be discarded.");
		outstream.println("touppercase=f      \t(tuc) Change all letters in reads and reference to upper-case.");
		outstream.println("                   \tValues: f (don't trim), r (trim right end), l (trim left end), n (convert to N instead of trimming).");
		outstream.println("qtrim=f            \tTrim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers. ");
		outstream.println("                   \tValues: t (trim both ends), f (neither end), r (right end only), l (left end only).");
		outstream.println("trimq=6            \tTrim quality threshold.");
		outstream.println("minlength=2        \t(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.");
		outstream.println("ziplevel=2         \t(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.");
		outstream.println("fastawrap=80       \tLength of lines in fasta output");
		outstream.println("qin=auto           \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto          \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
		outstream.println("rcomp=t            \tLook for reverse-complements of kmers also.");
	}
	
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Seal(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		FastaReadInputStream.SPLIT_READS=false;
		ByteFile.FORCE_MODE_BF2=Shared.threads()>2;
		
		/* Initialize local variables with defaults */
		boolean setOut=false, setOutb=false, qtrimRight_=false, qtrimLeft_=false;
		boolean rcomp_=true;
		boolean forbidNs_=false;
		boolean prealloc_=false;
		boolean useCountvector_=false;
		int tableType_=AbstractKmerTable.ARRAYH;
		int k_=31;
		int ways_=-1; //Currently disabled
		int minKmerHits_=1;
		long skipreads_=0;
		byte qin=-1, qout=-1;
		
		Parser parser=new Parser();
		parser.trimq=6;
		parser.minAvgQuality=0;
		parser.minReadLength=10;
		parser.maxReadLength=Integer.MAX_VALUE;
		parser.minLenFraction=0f;
		parser.requireBothBad=false;
		parser.maxNs=-1;
		boolean ordered_=false;
		int restrictLeft_=0, restrictRight_=0, speed_=0, qSkip_=1;
		int ambigMode_=AMBIG_RANDOM;
		int matchMode_=MATCH_ALL;
		boolean keepPairsTogether_=true;
		boolean printNonZeroOnly_=true;
		boolean rename_=false, useRefNames_=false;
		
		scaffoldNames.add(""); //Necessary so that the first real scaffold gets an id of 1, not zero
		scaffoldLengths.add(0);
		
		{
			boolean b=false;
			assert(b=true);
			EA=b;
		}
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseHist(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("out") || a.equals("out1") || a.equals("outm") || a.equals("outm1") || a.equals("outmatched") || a.equals("outmatched1")){
				outm1=b;
				setOut=true;
			}else if(a.equals("out2") || a.equals("outm") || a.equals("outm2") || a.equals("outmatched") || a.equals("outmatched2")){
				outm2=b;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("outunmatched") || a.equals("outunmatched1")){
				outu1=b;
			}else if(a.equals("outu2") || a.equals("outunmatched") || a.equals("outunmatched2")){
				outu2=b;
			}else if(a.equals("outs") || a.equals("outsingle")){
				outsingle=b;
			}else if(a.equals("stats") || a.equals("scafstats")){
				outstats=b;
			}else if(a.equals("refstats")){
				outrefstats=b;
			}else if(a.equals("rpkm") || a.equals("fpkm") || a.equals("cov") || a.equals("coverage")){
				outrpkm=b;
			}else if(a.equals("ref")){
				ref=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
			}else if(a.equals("literal")){
				literal=(b==null) ? null : b.split(",");
//				assert(false) : b+", "+Arrays.toString(literal);
			}else if(a.equals("forest")){
				boolean x=Tools.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.FOREST2D;}
			}else if(a.equals("array") || a.equals("array2")){
				boolean x=Tools.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAY2D;}
			}else if(a.equals("array1")){
				boolean x=Tools.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAY1D;}
			}else if(a.equals("arrayh") || a.equals("hybrid")){
				boolean x=Tools.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAYH;}
			}else if(a.equals("ways")){
				ways_=Integer.parseInt(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered_=Tools.parseBoolean(b);
				System.err.println("Set ORDERED to "+ordered_);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value greater than 0, such as k=27\n";
				k_=Integer.parseInt(b);
				assert(k_>0 && k_<32) : "k must be at least 1; default is 31.";
			}else if(a.equals("hdist") || a.equals("hammingdistance")){
				hammingDistance=Integer.parseInt(b);
				assert(hammingDistance>=0 && hammingDistance<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("edits") || a.equals("edist") || a.equals("editdistance")){
				editDistance=Integer.parseInt(b);
				assert(editDistance>=0 && editDistance<3) : "edit distance must be between 0 and 2; default is 0.";
			}else if(a.equals("skip") || a.equals("refskip") || a.equals("rskip")){
				refSkip=Integer.parseInt(b);
			}else if(a.equals("qskip")){
				qSkip_=Integer.parseInt(b);
			}else if(a.equals("speed")){
				speed_=Integer.parseInt(b);
				assert(speed_>=0 && speed_<=16) : "Speed range is 0 to 15.  Value: "+speed_;
			}else if(a.equals("skipreads")){
				skipreads_=Tools.parseKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("minkmerhits") || a.equals("minhits") || a.equals("mh")){
				minKmerHits_=Integer.parseInt(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				assert(WAYS>1) : "WAYS=1 is for debug mode.";
//				verbose=Tools.parseBoolean(b); //123
				if(verbose){outstream=System.err;} //For some reason System.out does not print in verbose mode.
			}else if(a.equals("mm") || a.equals("maskmiddle")){
				maskMiddle=Tools.parseBoolean(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("forbidns") || a.equals("forbidn") || a.equals("fn")){
				forbidNs_=Tools.parseBoolean(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isAlphabetic(b.charAt(0))){
					prealloc_=Tools.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc_=(preallocFraction>0);
				}
			}else if(a.equals("restrictleft")){
				restrictLeft_=Integer.parseInt(b);
			}else if(a.equals("restrictright")){
				restrictRight_=Integer.parseInt(b);
			}else if(a.equals("statscolumns") || a.equals("columns") || a.equals("cols")){
				STATS_COLUMNS=Integer.parseInt(b);
				assert(STATS_COLUMNS==3 || STATS_COLUMNS==5) : "statscolumns bust be either 3 or 5. Invalid value: "+STATS_COLUMNS;
			}else if(a.equals("ambiguous") || a.equals("ambig")){
				if(b==null){
					throw new RuntimeException(arg);
				}else if(b.equalsIgnoreCase("keep") || b.equalsIgnoreCase("best") || b.equalsIgnoreCase("first")){
					ambigMode_=AMBIG_FIRST;
				}else if(b.equalsIgnoreCase("all")){
					ambigMode_=AMBIG_ALL;
				}else if(b.equalsIgnoreCase("random") || b.equalsIgnoreCase("rand")){
					ambigMode_=AMBIG_RANDOM;
				}else if(b.equalsIgnoreCase("toss") || b.equalsIgnoreCase("discard") || b.equalsIgnoreCase("remove")){
					ambigMode_=AMBIG_TOSS;
				}else{
					throw new RuntimeException(arg);
				}
			}else if(a.equals("match") || a.equals("mode")){
				if(b==null){
					throw new RuntimeException(arg);
				}else if(b.equalsIgnoreCase("all") || b.equalsIgnoreCase("best")){
					matchMode_=MATCH_ALL;
				}else if(b.equalsIgnoreCase("first")){
					matchMode_=MATCH_FIRST;
				}else if(b.equalsIgnoreCase("unique") || b.equalsIgnoreCase("firstunique")){
					matchMode_=MATCH_UNIQUE;
				}else{
					throw new RuntimeException(arg);
				}
			}else if(a.equals("findbestmatch") || a.equals("fbm")){
				matchMode_=(Tools.parseBoolean(b) ? MATCH_ALL : MATCH_FIRST);
			}else if(a.equals("firstuniquematch") || a.equals("fum")){
				if(Tools.parseBoolean(b)){matchMode_=MATCH_UNIQUE;}
			}else if(a.equals("keeppairstogether") || a.equals("kpt")){
				keepPairsTogether_=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				printNonZeroOnly_=Tools.parseBoolean(b);
			}else if(a.equals("rename")){
				rename_=Tools.parseBoolean(b);
			}else if(a.equals("refnames") || a.equals("userefnames")){
				useRefNames_=Tools.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=(int)Tools.parseKMG(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("countvector")){
				useCountvector_=Tools.parseBoolean(b);
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && outu1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				outu1=args[i];
				setOut=true;
			}else if(i==2 && ref==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				ref=(new File(args[i]).exists() ? new String[] {args[i]} : args[i].split(","));
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Download parser fields
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
//			testsize=parser.testsize;
//			trimBadSequence=parser.trimBadSequence;
//			breakLength=parser.breakLength;
			
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			minLenFraction=parser.minLenFraction;
			minAvgQuality=parser.minAvgQuality;
			chastityFilter=parser.chastityFilter;
			minReadLength=parser.minReadLength;
			maxReadLength=parser.maxReadLength;
			maxNs=parser.maxNs;
//			untrim=parser.untrim;
//			minGC=parser.minGC;
//			maxGC=parser.maxGC;
//			filterGC=parser.filterGC;
//			minTrimLength=(parser.minTrimLength>=0 ? parser.minTrimLength : minTrimLength);
//			requireBothBad=parser.requireBothBad;
			removePairsIfEitherBad=!parser.requireBothBad;
		}
		
		if(ref!=null){
			for(String s : ref){refNames.add(s);}
		}
		if(literal!=null){refNames.add("literal");}
		refScafCounts=new int[refNames.size()];
		
		if(prealloc_){
			System.err.println("Note - if this program runs out of memory, please disable the prealloc flag.");
		}
		
		assert(outsingle==null) : "outsingle is currently not used.";
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		/* Set final variables; post-process and validate argument combinations */
		
		tableType=tableType_;
		hammingDistance=Tools.max(editDistance, hammingDistance);
		refSkip=Tools.max(0, refSkip);
		rcomp=rcomp_;
		forbidNs=(forbidNs_ || hammingDistance<1);
		skipreads=skipreads_;
		ORDERED=ordered_;
		restrictLeft=Tools.max(restrictLeft_, 0);
		restrictRight=Tools.max(restrictRight_, 0);
		ambigMode=ambigMode_;
		matchMode=matchMode_;
		keepPairsTogether=keepPairsTogether_;
		printNonZeroOnly=printNonZeroOnly_;
		rename=rename_;
		useRefNames=useRefNames_;
		speed=speed_;
		qSkip=qSkip_;
		noAccel=(speed<1 && qSkip<2);
		
		USE_COUNTVECTOR=useCountvector_;
		MAKE_QUALITY_HISTOGRAM=ReadStats.COLLECT_QUALITY_STATS;
		MAKE_QUALITY_ACCURACY=ReadStats.COLLECT_QUALITY_ACCURACY;
		MAKE_MATCH_HISTOGRAM=ReadStats.COLLECT_MATCH_STATS;
		MAKE_BASE_HISTOGRAM=ReadStats.COLLECT_BASE_STATS;
		MAKE_EHIST=ReadStats.COLLECT_ERROR_STATS;
		MAKE_INDELHIST=ReadStats.COLLECT_INDEL_STATS;
		MAKE_LHIST=ReadStats.COLLECT_LENGTH_STATS;
		MAKE_GCHIST=ReadStats.COLLECT_GC_STATS;
		MAKE_IDHIST=ReadStats.COLLECT_IDENTITY_STATS;
		
		if((speed>0 && qSkip>1) || (qSkip>1 && refSkip>1) || (speed>0 && refSkip>1)){
			System.err.println("WARNING: It is not recommended to use more than one of qskip, speed, and rskip together.");
			System.err.println("qskip="+qSkip+", speed="+speed+", rskip="+refSkip);
		}
		
		{
			long usableMemory;
			long tableMemory;

			{
				long memory=Runtime.getRuntime().maxMemory();
				double xmsRatio=Shared.xmsRatio();
				usableMemory=(long)Tools.max(((memory-96000000-(20*400000 /* for atomic arrays */))*(xmsRatio>0.97 ? 0.82 : 0.75)), memory*0.45);
				tableMemory=(long)(usableMemory*.95);
			}

			if(initialSize<1){
				final int factor=(tableType==AbstractKmerTable.ARRAY1D ? 12 : tableType==AbstractKmerTable.ARRAYH ? 22 : 27);
				final long memOverWays=tableMemory/(factor*WAYS);
				final double mem2=(prealloc_ ? preallocFraction : 1)*tableMemory;
				initialSize=(prealloc_ || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(factor*WAYS))) : initialSizeDefault);
				if(initialSize!=initialSizeDefault){
					System.err.println("Initial size set to "+initialSize);
				}
			}
		}
		
		k=k_;
		k2=k-1;
		minKmerHits=minKmerHits_;
		
		kfilter=(ref!=null || literal!=null);
		assert(kfilter==false || (k>0 && k<32)) : "K must range from 1 to 31.";
		
		middleMask=maskMiddle ? ~(3L<<(2*(k/2))) : -1L;
		
		
		/* Adjust I/O settings and filenames */
		
		if(qin!=-1 && qout!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY=false;
		}else if(qin!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.DETECT_QUALITY=false;
		}else if(qout!=-1){
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY_OUT=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(qfin1!=null && qfin1.contains("#") && in2!=null && !new File(qfin1).exists()){
			int pound=qfin1.lastIndexOf('#');
			String a=qfin1.substring(0, pound);
			String b=qfin1.substring(pound+1);
			qfin1=a+1+b;
			qfin2=a+2+b;
		}
		
		if(outu1!=null && outu1.contains("#")){
			int pound=outu1.lastIndexOf('#');
			String a=outu1.substring(0, pound);
			String b=outu1.substring(pound+1);
			outu1=a+1+b;
			outu2=a+2+b;
		}
		
		if(outm1!=null && outm1.contains("#")){
			int pound=outm1.lastIndexOf('#');
			String a=outm1.substring(0, pound);
			String b=outm1.substring(pound+1);
			outm1=a+1+b;
			outm2=a+2+b;
		}
		
		if((outu2!=null || outm2!=null) && (in1!=null && in2==null)){
			if(!FASTQ.FORCE_INTERLEAVED){System.err.println("Forcing interleaved input because paired output was specified for a single input file.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
		}

		if(!setOut){
			System.err.println("No output stream specified.  To write to stdout, please specify 'out=stdout.fq' or similar.");
//			out1="stdout";
//			outstream=System.err;
//			out2=null;
			outu1=outu2=null;
		}else if("stdout".equalsIgnoreCase(outu1) || "standarddout".equalsIgnoreCase(outu1)){
			outu1="stdout.fq";
			outstream=System.err;
			outu2=null;
		}

		if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, outm1, outm2, outsingle, outstats, outrpkm, outrefstats)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2, qfin1, qfin2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		if(!Tools.testInputFiles(true, true, ref)){
			throw new RuntimeException("\nCan't read to some reference files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, outu1, outu2, outm1, outm2, outsingle, outstats, outrpkm, outrefstats)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		assert(THREADS>0) : "THREADS must be greater than 0.";

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		
		if(ref==null && literal==null){
			System.err.println("ERROR: No reference sequences specified.  Use the -da flag to run anyway.");
			assert(false);
		}
				
		if(ref!=null){
			for(String s0 : ref){
				assert(s0!=null) : "Specified a null reference.";
				String s=s0.toLowerCase();
				assert(s==null || s.startsWith("stdin") || s.startsWith("standardin") || new File(s0).exists()) : "Can't find "+s0;
			}
		}
		
		//Initialize tables
		keySets=AbstractKmerTable.preallocate(WAYS, tableType, initialSize, (!prealloc_ || preallocFraction<1));
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void process(){
		
		/* Check for output file collisions */
		if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, outm1, outm2, outstats, outrpkm, outrefstats)){
			throw new RuntimeException("One or more output files were duplicate or could not be written to.  Check the names or set the 'overwrite=true' flag.");
		}
		
		/* Start overall timer */
		Timer t=new Timer();
		t.start();
		
//		boolean dq0=FASTQ.DETECT_QUALITY;
//		boolean ti0=FASTQ.TEST_INTERLEAVED;
//		int rbl0=Shared.READ_BUFFER_LENGTH;
//		FASTQ.DETECT_QUALITY=false;
//		FASTQ.TEST_INTERLEAVED=false;
//		Shared.READ_BUFFER_LENGTH=16;
		
		process2(t.time1);
		
//		FASTQ.DETECT_QUALITY=dq0;
//		FASTQ.TEST_INTERLEAVED=ti0;
//		Shared.READ_BUFFER_LENGTH=rbl0;
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		
		if(showSpeed){
			double rpnano=readsIn/(double)(t.elapsed);
			double bpnano=basesIn/(double)(t.elapsed);

			//Format with k or m suffixes
			String rpstring=(readsIn<100000 ? ""+readsIn : readsIn<100000000 ? (readsIn/1000)+"k" : (readsIn/1000000)+"m");
			String bpstring=(basesIn<100000 ? ""+basesIn : basesIn<100000000 ? (basesIn/1000)+"k" : (basesIn/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}

			outstream.println("\nTime:   \t\t\t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException("Seal terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(long startTime){
		
		/* Start phase timer */
		Timer t=new Timer();
		t.start();

		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			printMemory();
			outstream.println();
		}
		
		/* Fill tables with reference kmers */
		{
			final boolean oldTI=FASTQ.TEST_INTERLEAVED; //TODO: This needs to be changed to a non-static field, or somehow 'read mode' and 'ref mode' need to be distinguished.
			final boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			final boolean oldSplit=FastaReadInputStream.SPLIT_READS;
			final int oldML=FastaReadInputStream.MIN_READ_LEN;
			
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			FastaReadInputStream.SPLIT_READS=false;
			FastaReadInputStream.MIN_READ_LEN=1;
			
			storedKmers=spawnLoadThreads();
			
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			FastaReadInputStream.SPLIT_READS=oldSplit;
			FastaReadInputStream.MIN_READ_LEN=oldML;
			
			if(useRefNames){toRefNames();}
			t.stop();
		}
		
		/* Check memory */
		{
			long ram=freeMemory();
			ALLOW_LOCAL_ARRAYS=(scaffoldNames!=null && Tools.max(THREADS, 1)*3*8*scaffoldNames.size()<ram*5);
		}
		
		/* Dump kmers to text */
		if(dump!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(dump, overwrite, false, true);
			bsw.start();
			for(AbstractKmerTable set : keySets){
				set.dumpKmersAsBytes(bsw, k, 0);
			}
			bsw.poisonAndWait();
		}
		
		/* Do kmer matching of input reads */
		spawnProcessThreads(t);
		
		/* Unload kmers to save memory */
		if(RELEASE_TABLES){unloadKmers();}
		
		/* Write statistics to files */
		writeStats();
		writeRPKM();
		writeRefStats();
		
		/* Unload sequence data to save memory */
		if(RELEASE_TABLES){unloadScaffolds();}
		
		outstream.println("\nInput:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(ref!=null || literal!=null){
			outstream.println("Matched reads:          \t"+readsMatched+" reads ("+String.format("%.2f",readsMatched*100.0/readsIn)+"%) \t"+
					basesMatched+" bases ("+String.format("%.2f",basesMatched*100.0/basesIn)+"%)");
			outstream.println("Unmatched reads:        \t"+readsUnmatched+" reads ("+String.format("%.2f",readsUnmatched*100.0/readsIn)+"%) \t"+
					basesUnmatched+" bases ("+String.format("%.2f",basesUnmatched*100.0/basesIn)+"%)");
			outstream.flush();
		}
		if(qtrimLeft || qtrimRight){
			outstream.println("QTrimmed:               \t"+readsQTrimmed+" reads ("+String.format("%.2f",readsQTrimmed*100.0/readsIn)+"%) \t"+
					basesQTrimmed+" bases ("+String.format("%.2f",basesQTrimmed*100.0/basesIn)+"%)");
		}
		if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimModulo>0){
			outstream.println("FTrimmed:               \t"+readsFTrimmed+" reads ("+String.format("%.2f",readsFTrimmed*100.0/readsIn)+"%) \t"+
					basesFTrimmed+" bases ("+String.format("%.2f",basesFTrimmed*100.0/basesIn)+"%)");
		}
		if(minAvgQuality>0 || maxNs>=0){
			outstream.println("Low quality discards:   \t"+readsQFiltered+" reads ("+String.format("%.2f",readsQFiltered*100.0/readsIn)+"%) \t"+
					basesQFiltered+" bases ("+String.format("%.2f",basesQFiltered*100.0/basesIn)+"%)");
		}
//		outstream.println("Result:                 \t"+readsMatched+" reads ("+String.format("%.2f",readsMatched*100.0/readsIn)+"%) \t"+
//				basesMatched+" bases ("+String.format("%.2f",basesMatched*100.0/basesIn)+"%)");
	}
	
	/**
	 * Clear stored kmers.
	 */
	public void unloadKmers(){
		if(keySets!=null){
			for(int i=0; i<keySets.length; i++){keySets[i]=null;}
		}
	}
	
	/**
	 * Clear stored sequence data.
	 */
	public void unloadScaffolds(){
		if(scaffoldNames!=null && !scaffoldNames.isEmpty()){
			scaffoldNames.clear();
			scaffoldNames.trimToSize();
		}
		scaffoldReadCounts=null;
		scaffoldFragCounts=null;
		scaffoldBaseCounts=null;
		scaffoldLengths=null;
	}
	
	/**
	 * Write statistics about how many reads matched each reference scaffold.
	 */
	private void writeStats(){
		if(outstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outstats, overwrite, false, false);
		tsw.start();
		
		long rsum=0, bsum=0;
		
		/* Create StringNum list of scaffold names and hitcounts */
		ArrayList<StringNum> list=new ArrayList<StringNum>();
		for(int i=1; i<scaffoldNames.size(); i++){
			final long num1=scaffoldReadCounts.get(i), num2=scaffoldBaseCounts.get(i);
			if(num1>0 || !printNonZeroOnly){
				rsum+=num1;
				bsum+=num2;
				final String s=scaffoldNames.get(i);
				final int len=scaffoldLengths.get(i);
				final StringNum sn=new StringNum(s, len, num1, num2);
				list.add(sn);
			}
		}
		Collections.sort(list);
		final double rmult=100.0/(readsIn>0 ? readsIn : 1);
		final double bmult=100.0/(basesIn>0 ? basesIn : 1);
		
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		
		if(STATS_COLUMNS==3){
			tsw.print(String.format("#Total\t%d\n",readsIn));
			tsw.print(String.format("#Matched\t%d\t%.5f%%\n",rsum,rmult*rsum));
			tsw.print("#Name\tReads\tReadsPct\n");
			for(int i=0; i<list.size(); i++){
				StringNum sn=list.get(i);
				tsw.print(String.format("%s\t%d\t%.5f%%\n",sn.name,sn.reads,(sn.reads*rmult)));
			}
		}else{
			tsw.print(String.format("#Total\t%d\t%d\n",readsIn,basesIn));
			tsw.print(String.format("#Matched\t%d\t%.5f%%\n",rsum,rmult*rsum,bsum,bsum*bmult));
			tsw.print("#Name\tReads\tReadsPct\tBases\tBasesPct\n");
			for(int i=0; i<list.size(); i++){
				StringNum sn=list.get(i);
				tsw.print(String.format("%s\t%d\t%.5f%%\t%d\t%.5f%%\n",sn.name,sn.reads,(sn.reads*rmult),sn.bases,(sn.bases*bmult)));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write RPKM statistics.
	 */
	private void writeRPKM(){
		if(outrpkm==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrpkm, overwrite, false, false);
		tsw.start();

		/* Count mapped reads */
		long mapped=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mapped+=scaffoldReadCounts.get(i);
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(String.format("#Reads\t%d\n",readsIn));
		tsw.print(String.format("#Mapped\t%d\n",mapped));
		tsw.print(String.format("#RefSequences\t%d\n",Tools.max(0, scaffoldNames.size()-1)));
		tsw.print("#Name\tLength\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=1; i<scaffoldNames.size(); i++){
			final long reads=scaffoldReadCounts.get(i);
			final long frags=scaffoldFragCounts.get(i);
			final long bases=scaffoldBaseCounts.get(i);
			final String s=scaffoldNames.get(i);
			final int len=scaffoldLengths.get(i);
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(String.format("%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\n",s,len,bases,bases*invlen,reads,reads*mult2,frags,frags*mult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write statistics on a per-reference basis.
	 */
	private void writeRefStats(){
		if(outrefstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrefstats, overwrite, false, false);
		tsw.start();
		
		/* Count mapped reads */
		long mapped=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mapped+=scaffoldReadCounts.get(i);
		}
		
		final int numRefs=refNames.size();
		long[] refReadCounts=new long[numRefs];
		long[] refFragCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		long[] refLengths=new long[numRefs];
		
		for(int r=0, s=1; r<numRefs; r++){
			final int lim=s+refScafCounts[r];
			while(s<lim){
				refReadCounts[r]+=scaffoldReadCounts.get(s);
				refFragCounts[r]+=scaffoldFragCounts.get(s);
				refBaseCounts[r]+=scaffoldBaseCounts.get(s);
				refLengths[r]+=scaffoldLengths.get(s);
				s++;
			}
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(String.format("#Reads\t%d\n",readsIn));
		tsw.print(String.format("#Mapped\t%d\n",mapped));
		tsw.print(String.format("#References\t%d\n",Tools.max(0, refNames.size())));
		tsw.print("#Name\tLength\tScaffolds\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=0; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long frags=refFragCounts[i];
			final long bases=refBaseCounts[i];
			final long len=refLengths[i];
			final int scafs=refScafCounts[i];
			final String name=ReadWrite.stripToCore(refNames.get(i));
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(String.format("%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\n",name,len,scafs,bases,bases*invlen,reads,reads*mult2,frags,frags*mult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Fills the scaffold names array with reference names.
	 */
	private void toRefNames(){
		final int numRefs=refNames.size();
		for(int r=0, s=1; r<numRefs; r++){
			final int scafs=refScafCounts[r];
			final int lim=s+scafs;
			final String name=ReadWrite.stripToCore(refNames.get(r));
//			System.err.println("r="+r+", s="+s+", scafs="+scafs+", lim="+lim+", name="+name);
			while(s<lim){
//				System.err.println(r+", "+s+". Setting "+scaffoldNames.get(s)+" -> "+name);
				scaffoldNames.set(s, name);
				s++;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long spawnLoadThreads(){
		Timer t=new Timer();
		t.start();
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		/* For each reference file... */
		int refNum=0;
		if(ref!=null){
			for(String refname : ref){

				/* Start an input stream */
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, false, ff, null);
				cris.start(); //4567
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				/* Iterate through read lists from the input stream */
				while(reads!=null && reads.size()>0){
					{
						/* Assign a unique ID number to each scaffold */
						ArrayList<Read> reads2=new ArrayList<Read>(reads);
						for(Read r1 : reads2){
							final Read r2=r1.mate;
							final Integer id=scaffoldNames.size();
							refScafCounts[refNum]++;
							scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
							int len=r1.length();
							r1.obj=id;
							if(r2!=null){
								r2.obj=id;
								len+=r2.length();
							}
							scaffoldLengths.add(len);
						}

						/* Send a pointer to the read list to each LoadThread */
						for(LoadThread lt : loaders){
							boolean b=true;
							while(b){
								try {
									lt.queue.put(reads2);
									b=false;
								} catch (InterruptedException e) {
									//TODO:  This will hang due to still-running threads.
									throw new RuntimeException(e);
								}
							}
						}
					}

					/* Dispose of the old list and fetch a new one */
					cris.returnList(ln.id, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				/* Cleanup */
				cris.returnList(ln.id, ln.list.isEmpty());
				errorState|=ReadWrite.closeStream(cris);
				refNum++;
			}
		}

		/* If there are literal sequences to use as references */
		if(literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(literal.length);
			if(verbose){System.err.println("Adding literals "+Arrays.toString(literal));}

			/* Assign a unique ID number to each literal sequence */
			for(int i=0; i<literal.length; i++){
				final Integer id=scaffoldNames.size();
				final Read r=new Read(literal[i].getBytes(), null, id);
				refScafCounts[refNum]++;
				scaffoldNames.add(id.toString());
				scaffoldLengths.add(r.length());
				r.obj=id;
				list.add(r);
			}

			/* Send a pointer to the read list to each LoadThread */
			for(LoadThread lt : loaders){
				boolean b=true;
				while(b){
					try {
						lt.queue.put(list);
						b=false;
					} catch (InterruptedException e) {
						//TODO:  This will hang due to still-running threads.
						throw new RuntimeException(e);
					}
				}
			}
		}
		
		/* Signal loaders to terminate */
		for(LoadThread lt : loaders){
			boolean b=true;
			while(b){
				try {
					lt.queue.put(POISON);
					b=false;
				} catch (InterruptedException e) {
					//TODO:  This will hang due to still-running threads.
					throw new RuntimeException(e);
				}
			}
		}
		
		/* Wait for loaders to die, and gather statistics */
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					lt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			added+=lt.added;
			refKmers+=lt.refKmersT;
			refBases+=lt.refBasesT;
			refReads+=lt.refReadsT;
		}
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=WAYS;
		refBases/=WAYS;
		refReads/=WAYS;
		
		scaffoldReadCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldFragCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldBaseCounts=new AtomicLongArray(scaffoldNames.size());

		t.stop();
		if(DISPLAY_PROGRESS){
			outstream.println("Added "+added+" kmers; time: \t"+t);
			printMemory();
			outstream.println();
		}
		
		if(verbose){
			TextStreamWriter tsw=new TextStreamWriter("stdout", false, false, false, FileFormat.TEXT);
			tsw.start();
			for(AbstractKmerTable table : keySets){
				table.dumpKmersAsText(tsw, k, 1);
			}
			tsw.poisonAndWait();
		}
		
		return added;
	}
	
	/**
	 * Match reads against reference kmers, using multiple ProcessThread.
	 * @param t
	 */
	private void spawnProcessThreads(Timer t){
		t.start();
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, false, ff1, ff2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			cris.start(); //4567
		}
		final boolean paired=cris.paired();
		outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
		
		/* Create read output streams */
		final ConcurrentReadOutputStream rosm, rosu, ross;
		if(outu1!=null){
			final int buff=(!ORDERED ? 12 : Tools.max(32, 2*Shared.threads()));
			FileFormat ff1=FileFormat.testOutput(outu1, FileFormat.FASTQ, null, true, overwrite, append, ORDERED);
			FileFormat ff2=FileFormat.testOutput(outu2, FileFormat.FASTQ, null, true, overwrite, append, ORDERED);
			rosu=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, true);
			rosu.start();
		}else{rosu=null;}
		if(outm1!=null){
			final int buff=(!ORDERED ? 12 : Tools.max(32, 2*Shared.threads()));
			FileFormat ff1=FileFormat.testOutput(outm1, FileFormat.FASTQ, null, true, overwrite, append, ORDERED);
			FileFormat ff2=FileFormat.testOutput(outm2, FileFormat.FASTQ, null, true, overwrite, append, ORDERED);
			rosm=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, true);
			rosm.start();
		}else{rosm=null;}
		if(outsingle!=null){
			final int buff=(!ORDERED ? 12 : Tools.max(32, 2*Shared.threads()));
			FileFormat ff=FileFormat.testOutput(outsingle, FileFormat.FASTQ, null, true, overwrite, append, ORDERED);
			ross=ConcurrentReadOutputStream.getStream(ff, null, null, null, buff, null, true);
			ross.start();
		}else{ross=null;}
		if(rosu!=null || rosm!=null || ross!=null){
			t.stop();
			outstream.println("Started output streams:\t"+t);
			t.start();
		}
		
		/* Optionally skip the first reads, since initial reads may have lower quality */
		if(skipreads>0){
			long skipped=0;

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(skipped<skipreads && reads!=null && reads.size()>0){
				skipped+=reads.size();
				
				if(rosm!=null){rosm.add(new ArrayList<Read>(1), ln.id);}
				if(rosu!=null){rosu.add(new ArrayList<Read>(1), ln.id);}
				if(ross!=null){ross.add(new ArrayList<Read>(1), ln.id);}
				
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			if(reads==null || reads.isEmpty()){
				ReadWrite.closeStreams(cris, rosu, rosm, ross);
				System.err.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(cris, rosm, rosu, ross, ALLOW_LOCAL_ARRAYS));}
		for(ProcessThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ProcessThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			readsMatched+=pt.readsMatchedT;
			basesMatched+=pt.basesMatchedT;
			readsUnmatched+=pt.readsUnmatchedT;
			basesUnmatched+=pt.basesUnmatchedT;
			readsQTrimmed+=pt.readsQTrimmedT;
			basesQTrimmed+=pt.basesQTrimmedT;
			readsFTrimmed+=pt.readsFTrimmedT;
			basesFTrimmed+=pt.basesFTrimmedT;
			readsQFiltered+=pt.readsQFilteredT;
			basesQFiltered+=pt.basesQFilteredT;
			
			if(pt.scaffoldReadCountsT!=null && scaffoldReadCounts!=null){
				for(int i=0; i<pt.scaffoldReadCountsT.length; i++){scaffoldReadCounts.addAndGet(i, pt.scaffoldReadCountsT[i]);}
				pt.scaffoldReadCountsT=null;
			}
			if(pt.scaffoldBaseCountsT!=null && scaffoldBaseCounts!=null){
				for(int i=0; i<pt.scaffoldBaseCountsT.length; i++){scaffoldBaseCounts.addAndGet(i, pt.scaffoldBaseCountsT[i]);}
				pt.scaffoldBaseCountsT=null;
			}
			if(pt.scaffoldFragCountsT!=null && scaffoldFragCounts!=null){
				for(int i=0; i<pt.scaffoldFragCountsT.length; i++){scaffoldFragCounts.addAndGet(i, pt.scaffoldFragCountsT[i]);}
				pt.scaffoldFragCountsT=null;
			}
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris, rosu, rosm, ross);
		errorState|=ReadStats.writeAll(paired);
		
		t.stop();
		if(showSpeed){
			outstream.println("Processing time:   \t\t"+t);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	
	/**
	 * Loads kmers into a table.  Each thread handles all kmers X such that X%WAYS==tnum.
	 */
	private class LoadThread extends Thread{
		
		public LoadThread(final int tnum_){
			tnum=tnum_;
			map=keySets[tnum];
		}
		
		/**
		 * Get the next list of reads (or scaffolds) from the queue.
		 * @return List of reads
		 */
		private ArrayList<Read> fetch(){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			return list;
		}
		
		@Override
		public void run(){
			ArrayList<Read> reads=fetch();
			while(reads!=POISON){
				for(Read r1 : reads){
					assert(r1.pairnum()==0);
					final Read r2=r1.mate;

					final int rblen=(r1==null ? 0 : r1.length());
					final int rblen2=(r1.mateLength());
					
					added+=addToMap(r1, refSkip);
					if(r2!=null){
						added+=addToMap(r2, refSkip);
					}
				}
				reads=fetch();
			}
			
//			if(AbstractKmerTable.TESTMODE){
//				for(int i=0; i<ll.size; i++){
//					assert(map.contains(ll.get(i), il.get(i)));
//					assert(!map.contains(ll.get(i), Integer.MAX_VALUE));
//				}
//				ll=null;
//				il=null;
//			}
			
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
		}

		/**
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(final Read r, final int skip){
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			final long kmask=kMasks[k];
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;

			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return 0;}
			
			final int id=(Integer)r.obj;

			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=(rkmer>>>2)|(x2<<shift2);
					if(b=='N'){len=0;}else{len++;}
					if(verbose){System.err.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
							added+=addToMap(kmer, rkmer, k, extraBase, id, kmask);
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=(rkmer>>>2)|(x2<<shift2);
					if(b=='N'){len=0;}else{len++;}
					if(verbose){System.err.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
						final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask);
						added+=atm;
					}
				}
			}
			return added;
		}
		
		
		/**
		 * Adds this kmer to the table, including any mutations implied by editDistance or hammingDistance.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0
		 * @return Number of kmers stored
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0){
			
			assert(kmask0==kMasks[len]) : kmask0+", "+len+", "+kMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(kMasks[len]);
			
			if(verbose){System.err.println("addToMap_A; len="+len+"; kMasks[len]="+kMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hammingDistance==0){
				final long key=toValue(kmer, rkmer, kmask0);
				if(speed>0 && ((key/WAYS)&15)<speed){return 0;}
				if(key%WAYS!=tnum){return 0;}
				if(verbose){System.err.println("addToMap_B: "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+" = "+key);}
//				int[] old=map.getValues(key, new int[1]);
				added=map.set(key, id);
//				assert(old==null || map.contains(key, old));
//				assert(map.contains(key, id));
//				ll.add(key);
//				il.add(id); assert(AbstractKmerTable.TESTMODE);
				
//				if(AbstractKmerTable.TESTMODE){
//					for(int i=0; i<ll.size; i++){
//						assert(map.contains(ll.get(i), il.get(i)));
//						assert(!map.contains(ll.get(i), Integer.MAX_VALUE));
//					}
//				}
				
			}else if(editDistance>0){
//				long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, editDistance, extraBase);
			}else{
				added=mutate(kmer, rkmer, len, id, hammingDistance, -1);
			}
			if(verbose){System.err.println("addToMap added "+added+" keys.");}
			return added;
		}

//		private LongList ll=new LongList();
//		private IntList il=new IntList();
		
		/**
		 * Mutate and store this kmer through 'dist' recursions.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @param dist Number of mutations
		 * @param extraBase Base added to end in case of deletions
		 * @return Number of kmers stored
		 */
		private long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase){
			long added=0;
			
			final long key=toValue(kmer, rkmer, kMasks[len]);
			
			if(verbose){System.err.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+kMasks[len]);}
			if(key%WAYS==tnum){
				if(verbose){System.err.println("mutate_B: "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+" = "+key);}
				int x=map.set(key, id);
				if(verbose){System.err.println("mutate_B added "+x+" keys.");}
				added+=x;
				assert(map.contains(key));
			}
			
			if(dist>0){
				final int dist2=dist-1;
				
				//Sub
				for(int j=0; j<4; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
						if(temp!=kmer){
							long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
							added+=mutate(temp, rtemp, len, id, dist2, extraBase);
						}
					}
				}
				
				if(editDistance>0){
					//Del
					if(extraBase>=0 && extraBase<=3){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<2)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1);
							}
						}
					}

					//Ins
					final long eb2=kmer&3;
					for(int i=1; i<len; i++){
						final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>2);
						for(int j=0; j<4; j++){
							final long temp=temp0|setMasks[j][i-1];
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, eb2);
							}
						}
					}
				}

			}
			
			return added;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Number of kmers stored by this thread */
		public long added=0;
		/** Number of items encountered by this thread */
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		/** Thread number; used to determine which kmers to store */
		public final int tnum;
		/** Buffer of input read lists */
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);
		
		/** Destination for storing kmers */
		private final AbstractKmerTable map;
		
	}
	
	/*--------------------------------------------------------------*/

	/**
	 * Matches read kmers against reference kmers, performs binning and/or trimming, and writes output. 
	 */
	private class ProcessThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 * @param rosu_ Unmatched read output stream (optional)
		 * @param rosm_ Matched read output stream (optional)
		 */
		public ProcessThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosm_, ConcurrentReadOutputStream rosu_, ConcurrentReadOutputStream ross_, boolean localArrays){
			cris=cris_;
			rosm=rosm_;
			rosu=rosu_;
			ross=ross_;
			
			readstats=(MAKE_QUALITY_HISTOGRAM || MAKE_MATCH_HISTOGRAM || MAKE_BASE_HISTOGRAM || MAKE_EHIST || MAKE_INDELHIST || MAKE_LHIST || MAKE_GCHIST || MAKE_IDHIST) ? 
					new ReadStats() : null;
			
			final int alen=(scaffoldNames==null ? 0 : scaffoldNames.size());
			if(localArrays && alen>0 && alen<10000){
				scaffoldReadCountsT=new long[alen];
				scaffoldBaseCountsT=new long[alen];
				scaffoldFragCountsT=new long[alen];
			}else{
				scaffoldReadCountsT=scaffoldBaseCountsT=scaffoldFragCountsT=null;
			}
			
			if(USE_COUNTVECTOR){
				countVector=new IntList(1000);
				countArray=null;
			}else{
				countVector=null;
				countArray=new int[alen];
			}
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			ArrayList<Read> mlist=(rosm==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
			ArrayList<Read> ulist=(rosu==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
//			ArrayList<Read> single=(ross==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					if(readstats!=null){
						if(MAKE_QUALITY_HISTOGRAM){readstats.addToQualityHistogram(r1);}
						if(MAKE_BASE_HISTOGRAM){readstats.addToBaseHistogram(r1);}
						if(MAKE_MATCH_HISTOGRAM){readstats.addToMatchHistogram(r1);}
						if(MAKE_QUALITY_ACCURACY){readstats.addToQualityAccuracy(r1);}

						if(MAKE_EHIST){readstats.addToErrorHistogram(r1);}
						if(MAKE_INDELHIST){readstats.addToIndelHistogram(r1);}
						if(MAKE_LHIST){readstats.addToLengthHistogram(r1);}
						if(MAKE_GCHIST){readstats.addToGCHistogram(r1);}
						if(MAKE_IDHIST){readstats.addToIdentityHistogram(r1);}
					}

					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					readsInT++;
					basesInT+=initialLength1;
					if(r2!=null){
						readsInT++;
						basesInT+=initialLength2;
					}
					
					boolean remove=false;
					
					if(chastityFilter){
						if(r1!=null && r1.failsChastity()){
							basesQFilteredT+=r1.length();
							readsQFilteredT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && r2.failsChastity()){
							basesQFilteredT+=r2.length();
							readsQFilteredT++;
							r2.setDiscarded(true);
						}
					}
					
					if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimModulo>0){
						if(r1!=null && !r1.discarded()){
							final int len=r1.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b=Tools.min(b0, b1);
							final int x=TrimRead.trimToPosition(r1, a, b, 1);
							basesFTrimmedT+=x;
							readsFTrimmedT+=(x>0 ? 1 : 0);
							if(r1.length()<minlen1){r1.setDiscarded(true);}
						}
						if(r2!=null && !r2.discarded()){
							final int len=r2.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b=Tools.min(b0, b1);
							final int x=TrimRead.trimToPosition(r2, a, b, 1);
							basesFTrimmedT+=x;
							readsFTrimmedT+=(x>0 ? 1 : 0);
							if(r2.length()<minlen2){r2.setDiscarded(true);}
						}
					}
					
					if(removePairsIfEitherBad){remove=r1.discarded() || (r2!=null && r2.discarded());}
					else{remove=r1.discarded() && (r2==null || r2.discarded());}
					

					
					if(!remove){
						//Do quality trimming
						
						int rlen1=0, rlen2=0;
						if(r1!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
							}
							rlen1=r1.length();
							if(rlen1<minlen1 || rlen1>maxReadLength){
								r1.setDiscarded(true);
								if(verbose){System.err.println(r1.id+" discarded due to length.");}
							}
						}
						if(r2!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
							}
							rlen2=r2.length();
							if(rlen2<minlen2 || rlen2>maxReadLength){
								r2.setDiscarded(true);
								if(verbose){System.err.println(r2.id+" discarded due to length.");}
							}
						}
						
						//Discard reads if too short
						if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) || 
								(r1.discarded() && (r2==null || r2.discarded()))){
							basesQFilteredT+=(r1.length()+r1.mateLength());
							readsQTrimmedT+=1+r1.mateCount();
							remove=true;
						}
					}

					
					if(!remove){
						//Do quality filtering
						
						//Determine whether to discard the reads based on average quality
						if(minAvgQuality>0){
							if(r1!=null && r1.quality!=null && r1.avgQuality()<minAvgQuality){
								r1.setDiscarded(true);
								if(verbose){System.err.println(r1.id+" discarded due to low quality.");}
							}
							if(r2!=null && r2.quality!=null && r2.avgQuality()<minAvgQuality){
								r2.setDiscarded(true);
								if(verbose){System.err.println(r2.id+" discarded due to low quality.");}
							}
						}
						//Determine whether to discard the reads based on the presence of Ns
						if(maxNs>=0){
							if(r1!=null && r1.countUndefined()>maxNs){
								r1.setDiscarded(true);
								if(verbose){System.err.println(r1.id+" discarded due to Ns.");}
							}
							if(r2!=null && r2.countUndefined()>maxNs){
								r2.setDiscarded(true);
								if(verbose){System.err.println(r2.id+" discarded due to Ns.");}
							}
						}
						
						//Discard reads if too short
						if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) || 
								(r1.discarded() && (r2==null || r2.discarded()))){
							basesQFilteredT+=(r1.length()+r1.mateLength());
							readsQFilteredT+=1+r1.mateCount();
							remove=true;
						}
					}
					
					final int sites, assigned;
					if(remove){
						if(r1!=null){
							basesQFilteredT+=r1.length();
							readsQFilteredT++;
						}
						if(r2!=null){
							basesQFilteredT+=r2.length();
							readsQFilteredT++;
						}
						sites=assigned=0;
					}else{
						//Do kmer matching
						
						if(keepPairsTogether){
							
							final int a, b, max;
							if(countArray==null){
								countVector.size=0;
								a=findBestMatch(r1, keySets, countVector);
								b=findBestMatch(r2, keySets, countVector);
								if(verbose){System.err.println("countVector: "+countVector);}
								max=condenseLoose(countVector, idList1, countList1);
							}else{
								idList1.size=0;
								a=findBestMatch(r1, keySets, countArray, idList1);
								b=findBestMatch(r2, keySets, countArray, idList1);

								max=condenseLoose(countArray, idList1, countList1);
							}

							if(verbose){
								System.err.println("idList1: "+idList1);
								System.err.println("countList1: "+countList1);
							}
							if(rename){
								rename(r1, idList1, countList1);
								rename(r2, idList1, countList1);
							}
							filterTopScaffolds(idList1, countList1, finalList1, max);
							if(verbose){
								System.err.println("idList1: "+idList1);
								System.err.println("countList1: "+countList1);
								System.err.println("finalList1: "+finalList1);
							}
							sites=finalList1.size;
							
							if(max>=minKmerHits){
								assigned=assignTogether(r1, r2);
							}else{assigned=0;}
							
						}else{
							final int max1, max2, a, b;
							{
								if(countArray==null){
									countVector.size=0;
									a=findBestMatch(r1, keySets, countVector);
									max1=condenseLoose(countVector, idList1, countList1);
								}else{
									idList1.size=0;
									a=findBestMatch(r1, keySets, countArray, idList1);
									max1=condenseLoose(countArray, idList1, countList1);
								}
								if(rename){rename(r1, idList1, countList1);}
								filterTopScaffolds(idList1, countList1, finalList1, max1);
							}
							if(r2!=null){
								if(countArray==null){
									countVector.size=0;
									b=findBestMatch(r2, keySets, countVector);
									max2=condenseLoose(countVector, idList2, countList1);
								}else{
									idList2.size=0;
									b=findBestMatch(r2, keySets, countArray, idList2);
									max2=condenseLoose(countArray, idList2, countList2);
								}
								filterTopScaffolds(idList2, countList2, finalList2, max2);
								if(rename){rename(r2, idList2, countList2);}
							}else{max2=0;}
							
							sites=finalList1.size+finalList2.size;
							
							assigned=assignIndependently(r1, r2, max1, max2);
						}
					}
					
					if(remove || assigned<1){
						if(ulist!=null){ulist.add(r1);}
						
						//Track statistics
						if(r1!=null){
							readsUnmatchedT++;
							basesUnmatchedT+=r1.length();
						}
						if(r2!=null){
							readsUnmatchedT++;
							basesUnmatchedT+=r2.length();
						}
					}else{
						
						if(mlist!=null){mlist.add(r1);}
						
						//Track statistics
						if(r1!=null){
							readsMatchedT++;
							basesMatchedT+=r1.length();
						}
						if(r2!=null){
							readsMatchedT++;
							basesMatchedT+=r2.length();
						}
					}
				}
				
				//Send matched list to matched output stream
				if(rosu!=null){
					rosu.add(ulist, ln.id);
					ulist.clear();
				}
				
				//Send unmatched list to unmatched output stream
				if(rosm!=null){
					rosm.add(mlist, ln.id);
					mlist.clear();
				}
				
				if(ross!=null){
					assert(false);
//					ross.add(single, ln.id);
//					single.clear();
				}
				
				//Fetch a new read list
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}
		
		/**
		 * @param r
		 * @param idList
		 * @param countList
		 */
		private void rename(Read r, IntList idList, IntList countList) {
			if(r==null || idList.size<1){return;}
			StringBuilder sb=new StringBuilder();
			if(r.id==null){sb.append(r.numericID);}
			else{sb.append(r.id);}
			for(int i=0; i<idList.size; i++){
				int id=idList.get(i);
				int count=countList.get(i);
				sb.append('\t');
				sb.append(scaffoldNames.get(id));
				sb.append('=');
				sb.append(count);
			}
			r.id=sb.toString();
		}
		
		/**
		 * @param r1 Read 1
		 * @param r2 Read 2
		 * @return Number of sites assigned
		 */
		private int assignTogether(Read r1, Read r2){
			final int sites=finalList1.size;
			final int lenSum=r1.length()+(r1.mateLength());
			final int readSum=1+(r2==null ? 0 : 1);
			final int start, stop;

			if(sites<2 || ambigMode==AMBIG_ALL){
				start=0;
				stop=sites;
			}else if(ambigMode==AMBIG_TOSS){
				start=stop=0;
			}else if(ambigMode==AMBIG_FIRST){
				finalList1.sort();
				start=0;
				stop=1;
			}else if(ambigMode==AMBIG_RANDOM){
				start=(int)(r1.numericID%sites);
				stop=start+1;
			}else{
				throw new RuntimeException("Unknown mode "+ambigMode);
			}
			
			for(int j=start; j<stop; j++){
				int id=finalList1.get(j);
				
				if(scaffoldReadCountsT!=null){
					scaffoldReadCountsT[id]+=readSum;
					scaffoldBaseCountsT[id]+=lenSum;
					scaffoldFragCountsT[id]++;
				}else{
					scaffoldReadCounts.addAndGet(id, readSum);
					scaffoldBaseCounts.addAndGet(id, lenSum);
					scaffoldFragCounts.addAndGet(id, 1);
				}
			}
			return stop-start;
		}
		
		/**
		 * @param r1 Read 1
		 * @param r2 Read 2
		 * @param max1 Highest match count for read 1
		 * @param max2 Highest match count for read 2
		 * @return Number of sites assigned
		 */
		private int assignIndependently(Read r1, Read r2, int max1, int max2){
			int assigned=0;
			if(max1>=minKmerHits){
				final int sites=finalList1.size;
				final int lenSum=r1.length();
				final int start, stop;

				if(sites<2 || ambigMode==AMBIG_ALL){
					start=0;
					stop=sites;
				}else if(ambigMode==AMBIG_TOSS){
					start=stop=0;
				}else if(ambigMode==AMBIG_FIRST){
					finalList1.sort();
					start=0;
					stop=1;
				}else if(ambigMode==AMBIG_RANDOM){
					start=(int)(r1.numericID%sites);
					stop=start+1;
				}else{
					throw new RuntimeException("Unknown mode "+ambigMode);
				}
				
				for(int j=start; j<stop; j++){
					int id=finalList1.get(j);
					
					if(scaffoldReadCountsT!=null){
						scaffoldReadCountsT[id]++;
						scaffoldBaseCountsT[id]+=lenSum;
						if(max1>=max2){
							scaffoldFragCountsT[id]++;
						}
					}else{
						scaffoldReadCounts.addAndGet(id, 1);
						scaffoldBaseCounts.addAndGet(id, lenSum);
						if(max1>=max2){
							scaffoldFragCounts.addAndGet(id, 1);
						}
					}
				}
				assigned+=(stop-start);
			}
			
			if(max2>=minKmerHits){
				final int sites=finalList2.size;
				final int lenSum=r2.length();
				final int start, stop;

				if(sites<2 || ambigMode==AMBIG_ALL){
					start=0;
					stop=sites;
				}else if(ambigMode==AMBIG_TOSS){
					start=stop=0;
				}else if(ambigMode==AMBIG_FIRST){
					finalList2.sort();
					start=0;
					stop=1;
				}else if(ambigMode==AMBIG_RANDOM){
					start=(int)(r2.numericID%sites);
					stop=start+1;
				}else{
					throw new RuntimeException("Unknown mode "+ambigMode);
				}
				
				for(int j=start; j<stop; j++){
					int id=finalList2.get(j);
					
					if(scaffoldReadCountsT!=null){
						scaffoldReadCountsT[id]++;
						scaffoldBaseCountsT[id]+=lenSum;
						if(max2>max1){
							scaffoldFragCountsT[id]++;
						}
					}else{
						scaffoldReadCounts.addAndGet(id, 1);
						scaffoldBaseCounts.addAndGet(id, lenSum);
						if(max2>max1){
							scaffoldFragCounts.addAndGet(id, 1);
						}
					}
				}
				assigned+=(stop-start);
			}
			return assigned;
		}
		
		/**
		 * Pack a list of nonunique values into a list of unique values and a list of their counts.
		 * @param loose Nonunique values
		 * @param packed Unique values
		 * @param counts Counts of values
		 * @return
		 */
		private int condenseLoose(IntList loose, IntList packed, IntList counts){
			packed.size=0;
			counts.size=0;
			if(loose.size<1){return 0;}
			loose.sort();
			int prev=-1;
			int max=0;
			int count=0;
			for(int i=0; i<loose.size; i++){
				int id=loose.get(i);
//				System.err.println("i="+i+", id="+id+", count="+count+", prev="+prev);
				if(id==prev){
					count++;
				}else{
					if(count>0){
						packed.add(prev);
						counts.add(count);
						max=Tools.max(count, max);
//						assert(false) : "i="+i+", "+id+", "+count+", "+packed+", "+counts;
					}
					prev=id;
					count=1;
				}
			}
			if(count>0){
				packed.add(prev);
				counts.add(count);
				max=Tools.max(count, max);
			}
			return max;
		}
		
		/**
		 * Pack a list of counts from an array to an IntList.
		 * @param loose Counter array
		 * @param packed Unique values
		 * @param counts Counts of values
		 * @return
		 */
		private int condenseLoose(int[] loose, IntList packed, IntList counts){
			counts.size=0;
			if(packed.size<1){return 0;}

			int max=0;
			for(int i=0; i<packed.size; i++){
				final int p=packed.get(i);
				final int c=loose[p];
				counts.add(c);
				loose[p]=0;
				max=Tools.max(max, c);
			}
			return max;
		}
		
		private void filterTopScaffolds(IntList packed, IntList counts, IntList out, int thresh){
			out.size=0;
			if(packed.size<1){return;}
			for(int i=0; i<packed.size; i++){
				final int c=counts.get(i);
				assert(c>0 && c<=thresh) : c+"\n"+packed+"\n"+counts;
				if(c>=thresh){
					out.add(packed.get(i));
				}
			}
		}
		
		/**
		 * Returns number of matching kmers.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return number of total kmer matches
		 */
		private final int findBestMatch(final Read r, final AbstractKmerTable sets[], IntList hits){
			if(r==null || r.bases==null || storedKmers<1){return 0;}
			final byte[] bases=r.bases;
			final int minlen=k-1;
			final int minlen2=(maskMiddle ? k/2 : k);
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			final long kmask=kMasks[k];
			long kmer=0;
			long rkmer=0;
			int found=0;
			int len=0;
			
			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=Dedupe.baseToNumber[b];
				long x2=Dedupe.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(b=='N' && forbidNs){len=0;}else{len++;}
				if(verbose){System.err.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=minlen2 && i>=minlen){
					final long key=toValue(kmer, rkmer, kmask);
					if(noAccel || ((key/WAYS)&15)>=speed && (qSkip<1 || i%qSkip==0)){
						if(verbose){System.err.println("Testing key "+key);}
						AbstractKmerTable set=sets[(int)(key%WAYS)];
						if(verbose){System.err.println("Found set "+set.arrayLength()+", "+set.size());}
						final int[] ids=set.getValues(key, singleton);
						if(verbose){System.err.println("Found array "+(ids==null ? "null" : Arrays.toString(ids)));}
						if(ids!=null && ids[0]>-1){
							for(int id : ids){
								if(id==-1){break;}
								hits.add(id);
							}
							if(verbose){System.err.println("Found = "+(found+1)+"/"+minKmerHits);}
							found++;
//							assert(false) : (matchMode==MATCH_FIRST)+", "+(matchMode==MATCH_UNIQUE)+", "+ (ids.length==1 || ids[1]==-1);
							if(matchMode==MATCH_FIRST || (matchMode==MATCH_UNIQUE && (ids.length==1 || ids[1]==-1))){break;}
						}
					}
				}
			}
			return found;
		}
		
		/**
		 * Returns number of matching kmers.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return number of total kmer matches
		 */
		private final int findBestMatch(final Read r, final AbstractKmerTable sets[], int[] hits, IntList idList){
			if(r==null || r.bases==null || storedKmers<1){return 0;}
			final byte[] bases=r.bases;
			final int minlen=k-1;
			final int minlen2=(maskMiddle ? k/2 : k);
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			final long kmask=kMasks[k];
			long kmer=0;
			long rkmer=0;
			int found=0;
			int len=0;
			
			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=Dedupe.baseToNumber[b];
				long x2=Dedupe.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(b=='N' && forbidNs){len=0;}else{len++;}
				if(verbose){System.err.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=minlen2 && i>=minlen){
					final long key=toValue(kmer, rkmer, kmask);
					if(noAccel || ((key/WAYS)&15)>=speed && (qSkip<1 || i%qSkip==0)){
						if(verbose){System.err.println("Testing key "+key);}
						AbstractKmerTable set=sets[(int)(key%WAYS)];
						if(verbose){System.err.println("Found set "+set.arrayLength()+", "+set.size());}
						final int[] ids=set.getValues(key, singleton);
						if(verbose){System.err.println("Found array "+(ids==null ? "null" : Arrays.toString(ids)));}
						if(ids!=null && ids[0]>-1){
							for(int id : ids){
								if(id==-1){break;}
								hits[id]++;
								if(hits[id]==1){idList.add(id);}
							}
							if(verbose){System.err.println("Found = "+(found+1)+"/"+minKmerHits);}
							found++;
//							assert(false) : (matchMode==MATCH_FIRST)+", "+(matchMode==MATCH_UNIQUE)+", "+ (ids.length==1 || ids[1]==-1);
							if(matchMode==MATCH_FIRST || (matchMode==MATCH_UNIQUE && (ids.length==1 || ids[1]==-1))){break;}
						}
					}
				}
			}
			return found;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream cris;
		/** Output read streams */
		private final ConcurrentReadOutputStream rosm, rosu, ross;
		
		private final ReadStats readstats;
		private final IntList countVector;
		private final int[] countArray;
		
		private final IntList idList1=new IntList(), idList2=new IntList();
		private final IntList countList1=new IntList(), countList2=new IntList();
		private final IntList finalList1=new IntList(), finalList2=new IntList();
		
		long[] scaffoldReadCountsT;
		long[] scaffoldBaseCountsT;
		long[] scaffoldFragCountsT;
		final int[] singleton=new int[1];
		
		private long readsInT=0;
		private long basesInT=0;
		private long readsMatchedT=0;
		private long basesMatchedT=0;
		private long readsUnmatchedT=0;
		private long basesUnmatchedT=0;
		
		private long readsQTrimmedT=0;
		private long basesQTrimmedT=0;
		private long readsFTrimmedT=0;
		private long basesFTrimmedT=0;
		private long readsQFilteredT=0;
		private long basesQFilteredT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Object holding a String and numbers, for tracking the number of read and base hits per scaffold.
	 */
	private static class StringNum implements Comparable<StringNum>{
		
		public StringNum(String name_, int len_, long reads_, long bases_){
			name=name_;
			length=len_;
			reads=reads_;
			bases=bases_;
		}
		public final int compareTo(StringNum o){
			if(bases!=o.bases){return o.bases>bases ? 1 : -1;}
			if(reads!=o.reads){return o.reads>reads ? 1 : -1;}
			return name.compareTo(o.name);
		}
		public final boolean equals(StringNum o){
			return compareTo(o)==0;
		}
		public final String toString(){
			return name+"\t"+length+"\t"+reads+"\t"+bases;
		}
		
		/*--------------------------------------------------------------*/
		
		public final String name;
		public final int length;
		public final long reads, bases;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print statistics about current memory use and availability */
	private static final void printMemory(){
		if(GC_BEFORE_PRINT_MEMORY){
			System.gc();
			System.gc();
		}
		Runtime rt=Runtime.getRuntime();
		long mmemory=rt.maxMemory()/1000000;
		long tmemory=rt.totalMemory()/1000000;
		long fmemory=rt.freeMemory()/1000000;
		long umemory=tmemory-fmemory;
		outstream.println("Memory: "+/*"max="+mmemory+"m, total="+tmemory+"m, "+*/"free="+fmemory+"m, used="+umemory+"m");
	}
	
	/** Current available memory */
	private static final long freeMemory(){
		Runtime rt=Runtime.getRuntime();
		return rt.freeMemory();
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Fraction of available memory preallocated to arrays */
	private double preallocFraction=0.5;
	/** Initial size of data structures */
	private int initialSize=-1;
	/** Default initial size of data structures */
	private static final int initialSizeDefault=128000;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private final AbstractKmerTable[] keySets;
	/** A scaffold's name is stored at scaffoldNames.get(id).  
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	private final ArrayList<String> refNames=new ArrayList<String>();
	/** Number of scaffolds per reference. */
	private final int[] refScafCounts;
	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */ 
	private AtomicLongArray scaffoldReadCounts;
	/** scaffoldFragCounts[id] stores the number of fragments (reads or pairs) with kmer matches to that scaffold */ 
	private AtomicLongArray scaffoldFragCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */ 
	private AtomicLongArray scaffoldBaseCounts;
	/** Set to false to force threads to share atomic counter arrays. */ 
	private boolean ALLOW_LOCAL_ARRAYS=true;
	/** scaffoldLengths[id] stores the length of that scaffold */ 
	private IntList scaffoldLengths=new IntList();
	/** Array of reference files from which to load kmers */
	private String[] ref=null;
	/** Array of literal strings from which to load kmers */
	private String[] literal=null;
	
	/** Input reads */
	private String in1=null, in2=null;
	/** Input qual files */
	private String qfin1=null, qfin2=null;
	/** Output reads (matched and at least minlen) */
	private String outm1=null, outm2=null;
	/** Output reads (unmatched or shorter than minlen) */
	private String outu1=null, outu2=null;
	/** Output reads whose mate was unmatched or discarded */
	private String outsingle=null;
	/** Statistics output files */
	private String outstats=null, outrpkm=null, outrefstats=null;
	
	/** Dump kmers here. */
	private String dump=null;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	/** Process this fraction of input reads. */
	private float samplerate=1f;
	/** Set samplerate seed to this value. */
	private long sampleseed=-1;
	
	/** Output reads in input order.  May reduce speed. */
	private final boolean ORDERED;
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	private boolean maskMiddle=true;
	
	/** Store reference kmers with up to this many substitutions */
	private int hammingDistance=0;
	/** Store reference kmers with up to this many edits (including indels) */
	private int editDistance=0;
	/** Always skip this many kmers between used kmers when hashing reference. */
	private int refSkip=0;
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
	long readsIn=0;
	long basesIn=0;
	long readsMatched=0;
	long basesMatched=0;
	long readsUnmatched=0;
	long basesUnmatched=0;
	
	long readsQTrimmed=0;
	long basesQTrimmed=0;
	long readsFTrimmed=0;
	long basesFTrimmed=0;
	long readsQFiltered=0;
	long basesQFiltered=0;
	
	long refReads=0;
	long refBases=0;
	long refKmers=0;
	
	long storedKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Don't allow a read 'N' to match a reference 'A'.  
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	private final boolean forbidNs;
	/** AND bitmask with 0's at the middle base */ 
	private final long middleMask;
	/** Data structure to use.  Default: ARRAYH */
	private final int tableType;
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** A read must share at least this many kmers to be considered a match.  Default: 1 */
	private final int minKmerHits;
	/** Determines how to handle ambiguously-mapping reads */
	private final int ambigMode;
	/** Determines when to early-exit kmer matching */
	private final int matchMode;
	
	/** Quality-trim the left side */
	private final boolean qtrimLeft;
	/** Quality-trim the right side */
	private final boolean qtrimRight;
	/** Trim bases at this quality or below.  Default: 4 */
	private final byte trimq;
	/** Throw away reads below this average quality before trimming.  Default: 0 */
	private final byte minAvgQuality;
	/** Throw away reads failing chastity filter (:Y: in read header) */
	private final boolean chastityFilter;
	/** Throw away reads containing more than this many Ns.  Default: -1 (disabled) */
	private final int maxNs;
	/** Throw away reads shorter than this after trimming.  Default: 10 */
	private final int minReadLength;
	/** Throw away reads longer than this after trimming.  Default: Integer.MAX_VALUE */
	private final int maxReadLength;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	private final float minLenFraction;
	/** Filter reads by whether or not they have matching kmers */
	private final boolean kfilter;
	/** Trim left bases of the read to this position (exclusive, 0-based) */
	private final int forceTrimLeft;
	/** Trim right bases of the read after this position (exclusive, 0-based) */
	private final int forceTrimRight;
	/** Trim right bases of the read modulo this value. 
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	private final int forceTrimModulo;
	
	/** If positive, only look for kmer matches in the leftmost X bases */
	private int restrictLeft;
	/** If positive, only look for kmer matches the rightmost X bases */
	private int restrictRight;
	
	/** True iff java was launched with the -ea' flag */
	private final boolean EA;
	/** Skip this many initial input reads */
	private final long skipreads;

	/** Pairs go to outbad if either of them is bad, as opposed to requiring both to be bad.
	 * Default: true. */
	private final boolean removePairsIfEitherBad;
	
	/** Print only statistics for scaffolds that matched at least one read 
	 * Default: true. */ 
	private final boolean printNonZeroOnly;
	
	/** Rename reads to indicate what they matched.
	 * Default: false. */ 
	private final boolean rename;
	/** Use names of reference files instead of scaffolds.
	 * Default: false. */ 
	private final boolean useRefNames;
	
	/** Fraction of kmers to skip, 0 to 15 out of 16 */
	private final int speed;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int qSkip;
	
	/** True if speed and qSkip are disabled. */
	private final boolean noAccel;
	
	/** Pick a single scaffold per read pair, rather than per read */
	private final boolean keepPairsTogether;
	
	/** Store match IDs in an IntList rather than int array */ 
	private final boolean USE_COUNTVECTOR;
	
	private final boolean MAKE_QUALITY_ACCURACY;
	private final boolean MAKE_QUALITY_HISTOGRAM;
	private final boolean MAKE_MATCH_HISTOGRAM;
	private final boolean MAKE_BASE_HISTOGRAM;
	
	private final boolean MAKE_EHIST;
	private final boolean MAKE_INDELHIST;
	private final boolean MAKE_LHIST;
	private final boolean MAKE_GCHIST;
	private final boolean MAKE_IDHIST;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int VERSION=5;
	
	/** Number of tables (and threads, during loading) */ 
	private static final int WAYS=7; //123
	/** Verbose messages */
	public static final boolean verbose=false; //123
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	/** Indicates end of input stream */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;
	/** Number of columns for statistics output, 3 or 5 */
	public static int STATS_COLUMNS=5;
	/** Release memory used by kmer storage after processing reads */
	public static boolean RELEASE_TABLES=true;
	/** Max value of hitCount array */
	public static final int HITCOUNT_LEN=1000;
	
	/** x&clearMasks[i] will clear base i */
	private static final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	private static final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	private static final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	private static final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private static final long[] kMasks;
	
	public static HashMap<String,String> RQC_MAP=null;

	public static final int AMBIG_ALL=1, AMBIG_FIRST=2, AMBIG_TOSS=3, AMBIG_RANDOM=4;
	public static final int MATCH_ALL=1, MATCH_FIRST=2, MATCH_UNIQUE=3;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
	static{
		clearMasks=new long[32];
		leftMasks=new long[32];
		rightMasks=new long[32];
		kMasks=new long[32];
		setMasks=new long[4][32];
		for(int i=0; i<32; i++){
			clearMasks[i]=~(3L<<(2*i));
		}
		for(int i=0; i<32; i++){
			leftMasks[i]=((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			rightMasks[i]=~((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			kMasks[i]=((1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			for(long j=0; j<4; j++){
				setMasks[(int)j][i]=(j<<(2*i));
			}
		}
	}
	
}
