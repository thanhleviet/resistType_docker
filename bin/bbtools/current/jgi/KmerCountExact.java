package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.AbstractKmerTable;
import kmer.HashBuffer;
import kmer.KCountArray;
import kmer.KmerCount7MTA;
import kmer.KmerCountAbstract;
import kmer.Primes;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import dna.AminoAcid;
import dna.CoverageArray;
import dna.CoverageArray3;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import fileIO.FileFormat;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class KmerCountExact {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		KmerCountExact cke=new KmerCountExact(args);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		cke.process(t);
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx20g -cp <path> jgi.CountKmersExact in=<input file>");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("in2=<file>         \tUse this if 2nd read of pairs are in a different file.");
		outstream.println("out=<file>         \tDump kmers and counts to this file.");
		outstream.println("");
		outstream.println("threads=auto       \t(t) Set number of threads to use; default is number of logical processors.");
		outstream.println("overwrite=t        \t(ow) Set to false to force the program to abort rather than overwrite an existing file.");
		outstream.println("showspeed=t        \t(ss) Set to 'f' to suppress display of processing speed.");
		outstream.println("interleaved=auto   \t(int) If true, forces fastq input to be paired and interleaved.");
		outstream.println("k=28               \tKmer length used for finding contaminants.  Contaminants shorter than k will not be found.");
		outstream.println("minavgquality=0    \t(maq) Reads with average quality (before trimming) below this will be discarded.");
		outstream.println("touppercase=f      \t(tuc) Change all letters in reads and reference to upper-case.");
		outstream.println("qtrim=f            \tTrim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers. ");
		outstream.println("                   \tValues: t (trim both ends), f (neither end), r (right end only), l (left end only).");
		outstream.println("minq=4             \tTrim quality threshold.");
		outstream.println("minlength=2        \t(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.");
		outstream.println("ziplevel=2         \t(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.");
		outstream.println("fastawrap=80       \tLength of lines in fasta output");
		outstream.println("qin=auto           \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto          \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
		outstream.println("rcomp=t            \tLook for reverse-complements of kmers also.");
		outstream.println("forest=t           \tUse HashForest data structure");
		outstream.println("table=f            \tUse KmerTable data structure");
		outstream.println("array=f            \tUse HashArray data structure");
	}
	
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerCountExact(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		ByteFile.FORCE_MODE_BF2=Shared.threads()>2;
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		boolean setOut=false;
		boolean rcomp_=true;
		boolean useForest_=false, useTable_=false, useArray_=true, prealloc_=true;
		long skipreads_=0;
		int k_=31;
		int ways_=-1;
		int filterMax_=2;
		
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
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out") || a.equals("out1") || a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
				setOut=true;
			}else if(a.equals("mincounttodump") || a.equals("mindump") || a.equals("mincount")){
				minToDump=Integer.parseInt(b);
			}else if(a.equals("hist") || a.equals("khist")){
				outHist=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=Integer.parseInt(b);
			}else if(a.equals("forest")){
				useForest_=Tools.parseBoolean(b);
				if(useForest_){useTable_=useArray_=false;}
			}else if(a.equals("table")){
				useTable_=Tools.parseBoolean(b);
				if(useTable_){useForest_=useArray_=false;}
			}else if(a.equals("array")){
				useArray_=Tools.parseBoolean(b);
				if(useArray_){useTable_=useForest_=false;}
			}else if(a.equals("ways")){
				ways_=Integer.parseInt(b);
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
				buflen=Integer.parseInt(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value from 1 to 31, such as k=28\n";
				k_=Integer.parseInt(b);
			}else if(a.equals("skipreads")){
				skipreads_=Tools.parseKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				verbose=Tools.parseBoolean(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isAlphabetic(b.charAt(0))){
					prealloc_=Tools.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc_=(preallocFraction>0);
				}
			}else if(a.equals("prefilter")){
				if(b==null || b.length()<1 || !Character.isDigit(b.charAt(0))){
					prefilter=Tools.parseBoolean(b);
				}else{
					filterMax_=Integer.parseInt(b);
					prefilter=filterMax_>0;
				}
			}else if(a.equals("prefiltersize")){
				prefilterFraction=Tools.max(0, Double.parseDouble(b));
				assert(prefilterFraction<=1) : "prefiltersize must be 0-1, a fraction of total memory.";
				prefilter=prefilterFraction>0;
			}else if(a.equals("prehashes") || a.equals("hashes")){
				prehashes=Integer.parseInt(b);
			}else if(a.equals("onepass")){
				onePass=Tools.parseBoolean(b);
			}else if(a.equals("passes")){
				int passes=Integer.parseInt(b);
				onePass=(passes<2);
			}else if(a.equals("histcolumns")){
				histColumns=Integer.parseInt(b);
			}else if(a.equals("histmax")){
				histMax=Integer.parseInt(b);
			}else if(a.equals("histheader")){
				histHeader=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				histZeros=!Tools.parseBoolean(b);
			}else if(a.equals("minheight")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak")){
				maxPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(a.equals("peaks") || a.equals("peaksout")){
				outPeaks=b;
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && outKmers==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				outKmers=args[i];
				setOut=true;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Download parser fields
//			forceTrimLeft=parser.forceTrimLeft;
//			forceTrimRight=parser.forceTrimRight;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
//			minLenFraction=parser.minLenFraction;
			minAvgQuality=parser.minAvgQuality;
//			minReadLength=parser.minReadLength;
//			maxReadLength=parser.maxReadLength;
//			maxNs=parser.maxNs;
//			untrim=parser.untrim;
//			minGC=parser.minGC;
//			maxGC=parser.maxGC;
//			filterGC=parser.filterGC;
//			minTrimLength=(parser.minTrimLength>=0 ? parser.minTrimLength : minTrimLength);
//			requireBothBad=parser.requireBothBad;
//			removePairsIfEitherBad=!parser.requireBothBad;
		}
		
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
//			long tmemory=Runtime.getRuntime().totalMemory();
			usableMemory=(long)Tools.max(((memory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.75)), memory*0.45);
			filterMemory=(long)(prefilter ? usableMemory*prefilterFraction : 0);
			tableMemory=(long)(usableMemory*.95-filterMemory);
		}
		
		if(ways_<1){
			long maxKmers=(2*tableMemory)/12;
			long minWays=Tools.min(10000, maxKmers/Integer.MAX_VALUE);
			ways_=(int)Tools.max(31, THREADS*2, minWays);
			ways_=(int)Primes.primeAtLeast(ways_);
			assert(ways_>0);
			System.err.println("ways="+ways_);
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		onePass=onePass&prefilter;
		prealloc=prealloc_;
		useForest=useForest_;
		useTable=useTable_;
		useArray=useArray_;
		rcomp=rcomp_;
		skipreads=skipreads_;
		WAYS=ways_;
		filterMax=Tools.min(filterMax_, 0x7FFFFFFF);
		
		k=k_;
		k2=k-1;
		
		if(initialSize<1){
			final long memOverWays=tableMemory/(12*WAYS);
			final double mem2=(prealloc ? preallocFraction : 1)*tableMemory;
			initialSize=(prealloc || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(12*WAYS))) : initialSizeDefault);
			if(initialSize!=initialSizeDefault){
				System.err.println("Initial size set to "+initialSize);
			}
		}
		
		/* Adjust I/O settings and filenames */
		
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
		
		if(outKmers!=null && !Tools.canWrite(outKmers, overwrite)){throw new RuntimeException("Output file "+outKmers+" already exists, and overwrite="+overwrite);}

		assert(!in1.equalsIgnoreCase(outKmers));
		assert(!in1.equalsIgnoreCase(in2));
		assert(THREADS>0);

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		
		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			printMemory();
			outstream.println();
		}
		
		final int tableType=(useForest ? AbstractKmerTable.FOREST1D : useTable ? AbstractKmerTable.TABLE : useArray ? AbstractKmerTable.ARRAY1D : 0);
		keySets=AbstractKmerTable.preallocate(WAYS, tableType, initialSize, (!prealloc || preallocFraction<1));
		
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void process(Timer t){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, outKmers, outHist, outPeaks);
		
		/* Count kmers */
		process2();
		
		if(THREADS>1 && (outHist!=null || outPeaks!=null) && outKmers!=null){
			Timer tout=new Timer();
			tout.start();
			Thread a=new DumpKmersThread();
			Thread b=new MakeKhistThread();
			a.start();
			b.start();
			while(a.getState()!=Thread.State.TERMINATED){
				try {
					a.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			while(b.getState()!=Thread.State.TERMINATED){
				try {
					b.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			tout.stop();
			outstream.println("Write Time:                 \t"+tout);
		}else{
			if(outHist!=null){
				makeKhist(outHist, outPeaks, histColumns, histMax, histHeader, histZeros, true);
			}
			if(outKmers!=null){
				dumpKmersAsText(outKmers, k, minToDump, true);
			}
		}
		
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
			
			outstream.println("Total Time:               \t"+t);
			outstream.println("\nReads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException("BBDuk terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(){
		
		/* Start phase timer */
		Timer t=new Timer();
		t.start();

		if(DISPLAY_PROGRESS){
			outstream.println("Before Load:");
			printMemory();
			outstream.println();
		}
		
		prefilterArray=null;
//		outstream.println();
		if(prefilter){
			KmerCountAbstract.CANONICAL=true;
			
			long precells=-1;
			int cbits=1;
			if(onePass){
				while(filterMax>=(1<<cbits)){cbits*=2;}
			}else{
				while(filterMax+1>=(1<<cbits)){cbits*=2;}
			}
			byte minq=0;
			if(precells<1){
				long prebits=(filterMemory-10)*8;
				precells=prebits/cbits;
			}
			if(prehashes<1){prehashes=2;}
			
			if(onePass){
				prefilterArray=KmerCount7MTA.makeKca(null, null, null, k, cbits, 0, precells, prehashes, minq, true, maxReads, 1, 1, 1, 1, null);
			}else if(precells>100000){
				Timer ht=new Timer();
				ht.start();
				
				ArrayList<String> extra=null;
				prefilterArray=KmerCount7MTA.makeKca(in1, in2, extra, k, cbits, 0, precells, prehashes, minq, true, maxReads, 1, 1, 1, 1, null);
				assert(filterMax<prefilterArray.maxValue);
				outstream.println("Made prefilter:   \t"+prefilterArray.toShortString(prehashes));
				double uf=prefilterArray.usedFraction();
				if(uf>0.6){
					outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" : 
						uf>0.7 ? "fairly" : "somewhat")+" full.  Ideal load is under 60% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; or increase prefiltersize.");
				}
				ht.stop();

				if(DISPLAY_PROGRESS){
					outstream.println("Prefilter time:\t"+ht);
					outstream.println("After prefilter:");
					printMemory();
					outstream.println();
				}
			}else{prefilter=false;}
		}
		
		/* Fill tables with kmers */
		long added=loadKmers();
		
		if(DISPLAY_PROGRESS){
			outstream.println("Final:");
			printMemory();
			outstream.println();
		}
		
		t.stop();
		
		/* Write statistics to files */
//		writeStats(System.nanoTime()-startTime);
		
		outstream.println("Input:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(qtrimLeft || qtrimRight){
			outstream.println("QTrimmed:               \t"+readsTrimmed+" reads ("+String.format("%.2f",readsTrimmed*100.0/readsIn)+"%) \t"+
					basesTrimmed+" bases ("+String.format("%.2f",basesTrimmed*100.0/basesIn)+"%)");
		}
		if(minAvgQuality>0){
			outstream.println("Low quality discards:   \t"+lowqReads+" reads ("+String.format("%.2f",lowqReads*100.0/readsIn)+"%) \t"+
					lowqBases+" bases ("+String.format("%.2f",lowqBases*100.0/basesIn)+"%)");
		}
//		outstream.println("Result:                 \t"+readsOut+" reads ("+String.format("%.2f",readsOut*100.0/readsIn)+"%) \t"+
//				basesOut+" bases ("+String.format("%.2f",basesOut*100.0/basesIn)+"%)");
		outstream.println("Unique Kmers:               \t"+added);
		
		outstream.println("\nTable Load Time:          \t"+t);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	/**
	 * load reads into tables, using multiple ProcessThread.
	 */
	private long loadKmers(){
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, false, ff1, ff2);
			cris.start(); //4567
		}
		
		/* Optionally skip the first reads, since initial reads may have lower quality */
		if(skipreads>0){
			long skipped=0;

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(skipped<skipreads && reads!=null && reads.size()>0){
				skipped+=reads.size();
				
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			if(reads==null || reads.isEmpty()){
				ReadWrite.closeStreams(cris);
				System.err.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(cris));}
		for(ProcessThread pt : alpt){pt.start();}
		
		long added=0;
		
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
			added+=pt.added;
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			readsTrimmed+=pt.readsTrimmedT;
			basesTrimmed+=pt.basesTrimmedT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Matches read kmers against reference kmers, performs binning and/or trimming, and writes output. 
	 */
	private class ProcessThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 * @param ros_ Unmatched read output stream (optional)
		 * @param rosb_ Matched read output stream (optional)
		 */
		public ProcessThread(ConcurrentReadInputStream cris_){
			cris=cris_;
			table=new HashBuffer(keySets, buflen);
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					readsInT++;
					basesInT+=r1.length();
					if(r2!=null){
						readsInT++;
						basesInT+=r2.length();
					}
					
					//Determine whether to discard the reads based on average quality
					if(minAvgQuality>0){
						if(r1!=null && r1.quality!=null && r1.avgQuality()<minAvgQuality){r1.setDiscarded(true);}
						if(r2!=null && r2.quality!=null && r2.avgQuality()<minAvgQuality){r2.setDiscarded(true);}
					}
					
					int rlen1=0, rlen2=0;
					if(r1!=null){
						if(qtrimLeft || qtrimRight){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						rlen1=r1.length();
						if(rlen1<k){r1.setDiscarded(true);}
					}
					if(r2!=null){
						if(qtrimLeft || qtrimRight){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						rlen2=r2.length();
						if(rlen2<k){r2.setDiscarded(true);}
					}

					if(r1!=null){
						if(r1.discarded()){
							lowqBasesT+=r1.length();
							lowqReadsT++;
						}else{
							long temp=addKmersToTable(r1);
							added+=temp;
							if(verbose){System.err.println("Added "+temp);}
						}
					}
					if(r2!=null){
						if(r2.discarded()){
							lowqBasesT+=r2.length();
							lowqReadsT++;
						}else{
							long temp=addKmersToTable(r2);
							added+=temp;
							if(verbose){System.err.println("Added "+temp);}
						}
					}
				}
				
				//Fetch a new read list
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			added+=table.flush();
		}
		
		

		
		/**
		 * Counts the number of kmer hits for a read.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return Number of hits
		 */
		private final int addKmersToTable(final Read r){
			if(onePass){return addKmersToTable_onePass(r);}
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final long key=toValue(kmer, rkmer);
					if(!prefilter || prefilterArray.read(key)>filterMax){
						int temp=table.incrementAndReturnNumCreated(key);
						created+=temp;
						if(verbose){System.err.println("Added "+temp);}
					}
				}
			}
			return created;
		}
		
		

		
		/**
		 * Counts the number of kmer hits for a read.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return Number of hits
		 */
		private final int addKmersToTable_onePass(final Read r){
			assert(prefilter);
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final long key=toValue(kmer, rkmer);
					int count=prefilterArray.incrementAndReturnUnincremented(key, 1);
					if(count>=filterMax){
						int temp=table.incrementAndReturnNumCreated(key);
						created+=temp;
						if(verbose){System.err.println("Added "+temp);}
					}
				}
			}
			return created;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream cris;
		
		private final HashBuffer table;
		
		public long added=0;
		
		private long readsInT=0;
		private long basesInT=0;
		private long readsOutT=0;
		private long basesOutT=0;
		private long readsTrimmedT=0;
		private long basesTrimmedT=0;
		private long lowqReadsT=0;
		private long lowqBasesT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public boolean dumpKmersAsText(String fname, int k, int minToDump, boolean printTime){
		if(fname==null){return false;}
		Timer t=new Timer();
		t.start();
		
//		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, true);
//		tsw.start();
//		for(AbstractKmerTable set : keySets){
//			set.dumpKmersAsText(tsw, k, minToDump);
//		}
//		tsw.poisonAndWait();
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();
		for(AbstractKmerTable set : keySets){
			set.dumpKmersAsBytes(bsw, k, minToDump);
		}
		bsw.poisonAndWait();
		
		t.stop();
		if(printTime){outstream.println("Kmer Dump Time:             \t"+t);}
		return bsw.errorState;
	}
	
	public boolean makeKhist(String fname, String peaks, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime){
		if(fname==null){return false;}
		Timer t=new Timer();
		t.start();
		
//		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, true);
//		tsw.start();
//		if(printHeader){
//			tsw.print("#Depth\t"+(cols==3 ? "RawCount\t" : "")+"Count\n");
//		}
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();
		if(printHeader){
			bsw.print("#Depth\t"+(cols==3 ? "RawCount\t" : "")+"Count\n");
		}
		
		CoverageArray ca=new CoverageArray3();
		for(AbstractKmerTable set : keySets){
			set.fillHistogram(ca, histMax);
		}
		
		if(peaks!=null){
			CallPeaks.printClass=false;
			long[] array=Tools.toArray(ca);
			CallPeaks.printPeaks(array, peaks, overwrite, minHeight, minVolume, minWidth, minPeak, maxPeak, maxPeakCount);
		}
		
//		StringBuilder sb=new StringBuilder();
		for(int i=1; i<=ca.maxIndex; i++){
			int count=ca.get(i);
			if(printZeros || count>0){
//				sb.append(i);
//				sb.append('\t');
//				if(cols==3){
//					sb.append(i*(long)count);
//					sb.append('\t');
//				}
//				sb.append(count);
//				sb.append('\n');
//				bsw.print(sb.toString());
//				sb.setLength(0);

				bsw.print(i);
				bsw.print('\t');
				if(cols==3){
					bsw.print(i*(long)count);
					bsw.print('\t');
				}
				bsw.print(count);
				bsw.print('\n');
			}
		}
		bsw.poisonAndWait();
		t.stop();
		if(printTime){outstream.println("Histogram Write Time:       \t"+t);}
		return bsw.errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class DumpKmersThread extends Thread {
		
		DumpKmersThread(){}
		
		public void run(){
			dumpKmersAsText(outKmers, k, minToDump, false);
		}
		
	}
	
	private class MakeKhistThread extends Thread {
		
		MakeKhistThread(){}
		
		public void run(){
			makeKhist(outHist, outPeaks, histColumns, histMax, histHeader, histZeros, false);
		}
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
		outstream.println("Memory: "+"max="+mmemory+"m, total="+tmemory+"m, "+"free="+fmemory+"m, used="+umemory+"m");
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	private final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Has this class encountered errors while processing? */
	public boolean errorState=false;

	/** Use a count-min prefilter for low-depth kmers */
	public boolean prefilter=false;
	/** Fill the prefilter at the same time as the main table */
	public boolean onePass=false;
	/** Number of hashes used by prefilter */
	public int prehashes=2;
	/** Fraction of memory used by prefilter */
	private double prefilterFraction=0.2;
	
	/** Initial size of data structures */
	private int initialSize=-1;
	private static final int initialSizeDefault=128000;
	/** Fraction of available memory preallocated to arrays */
	private double preallocFraction=1.0;
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private final AbstractKmerTable[] keySets;
	/** A scaffold's name is stored at scaffoldNames.get(id).  
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private final ArrayList<String> scaffoldNames=new ArrayList<String>();
	
	private KCountArray prefilterArray=null;
	
	/** Input reads */
	private String in1=null, in2=null;
	/** Kmer count output file */
	private String outKmers=null;
	/** Histogram output file */
	private String outHist=null;
	/** Histogram peak output file */
	private String outPeaks=null;

	/** Histogram columns */
	private int histColumns=2;
	/** Histogram rows */
	private int histMax=100000;
	/** Histogram columns */
	private boolean histHeader=false;
	/** Histogram show rows with 0 count */
	private boolean histZeros=false;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	
	private int buflen=1000;
	
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	long readsTrimmed=0;
	long basesTrimmed=0;
	long lowqReads=0;
	long lowqBases=0;

	private long minHeight=2;
	private long minVolume=2;
	private int minWidth=2;
	private int minPeak=2;
	private int maxPeak=Integer.MAX_VALUE;
	private int maxPeakCount=8;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	private final long usableMemory;
	private final long filterMemory;
	private final long tableMemory;
	
	/** Number of tables (and threads, during loading) */ 
	private final boolean prealloc;
	
	/** Number of tables (and threads, during loading) */ 
	private final int WAYS;
	
	/** Filter kmers up to this level; don't store them in primary data structure */ 
	private final int filterMax;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Use HashForest data structure */
	private final boolean useForest;
	/** Use KmerTable data structure */
	private final boolean useTable;
	/** Use HashArray data structure (default) */
	private final boolean useArray;	
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	
	/** min kmer count to dump to text */
	private int minToDump=1;
	
	/** Quality-trim the left side */
	private final boolean qtrimLeft;
	/** Quality-trim the right side */
	private final boolean qtrimRight;
	/** Trim bases at this quality or below.  Default: 4 */
	private final byte trimq;
	/** Throw away reads below this average quality before trimming.  Default: 0 */
	private final byte minAvgQuality;
	
	/** True iff java was launched with the -ea' flag */
	private final boolean EA;
	/** Skip this many initial input reads */
	private final long skipreads;
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int VERSION=1;
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=false;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
//	static{
//	}
	

	
}
