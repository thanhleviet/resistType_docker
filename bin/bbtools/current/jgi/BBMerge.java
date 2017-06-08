package jgi;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.KCountArray;
import kmer.KmerCount7MTA;
import kmer.KmerCountAbstract;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;

import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;

import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;

/**
 * @author Brian Bushnell
 * @date Aug 14, 2012
 *
 */
public class BBMerge {
	
	
	public static void main(String[] args){
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		BBMerge mr=new BBMerge(args);
		mr.process();
	}
	
	
	private static void printOptions(){
		System.err.println("Syntax:\n");
		System.err.println("java -ea -Xmx200m -cp <path> jgi.BBMerge in=<input file> out=<output file>");
		System.err.println("For more options, please run the shellscript with no arguments, or look at its contents.");
	}
	
	
	private static String[] preparse(String[] args){
		if(args==null){return new String[0];}
		int nulls=0;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);}
			if(a.equals("jni") || a.equals("usejni")){
				Shared.USE_JNI=Tools.parseBoolean(b);
			}else if(a.equals("strict")){
				strict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("loose")){
				loose=Tools.parseBoolean(b);
				if(loose){strict=false;}
				args[i]=null;
				nulls++;
			}else if(a.equals("vloose")){
				vloose=Tools.parseBoolean(b);
				if(vloose){strict=false;}
				args[i]=null;
				nulls++;
			}else if(a.equals("fast")){
				fast=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("normal")){
				if(Tools.parseBoolean(b)){
					strict=loose=vloose=fast=false;
				}
				args[i]=null;
				nulls++;
			}
		}
		if(nulls==0){return args;}
		ArrayList<String> args2=new ArrayList<String>(args.length-nulls+5);
		if(strict){
			loose=vloose=false;
			args2.add("minoverlap0=3");
			args2.add("maxbad=3");
			args2.add("margin=3");
			args2.add("minqo=10");
			args2.add("qualiters=2");
			args2.add("minoverlap=14");
			args2.add("minentropy=56");
		}else if(vloose){
			loose=true;
			strict=false;
			args2.add("minoverlap=8");
			args2.add("minoverlap0=9");
			args2.add("minoi=35");
			args2.add("qualiters=4");
			args2.add("minentropy=26");
			args2.add("mismatches=4");
			args2.add("margin=2");
			args2.add("efilter=t");
		}else if(loose){
			strict=false;
			args2.add("minoverlap=8");
			args2.add("minoverlap0=9");
			args2.add("minoi=35");
			args2.add("qualiters=4");
			args2.add("minentropy=30");
			args2.add("mismatches=3");
			args2.add("margin=2");
			args2.add("efilter=t");
		}
		
		if(fast){
			args2.add("qualiters=1");
			args2.add("tof=0");
		}
		
		for(String s : args){
			if(s!=null){args2.add(s);}
		}
		return args2.toArray(new String[args2.size()]);
	}
	
	public BBMerge(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		System.err.println("BBMerge version "+version);
		
		args=preparse(args);
		
		Timer ttotal=new Timer();
		ttotal.start();
		
		in1=(args[0].indexOf('=')>0 ? null : args[0]);
		in2=(in1!=null && args.length>1 && args[1].indexOf('=')<0 ? args[1] : null);
		if(in2!=null && "null".equalsIgnoreCase(in2)){in2=null;}
		
		{
			if(in1!=null && !in1.contains(",") && !in1.startsWith("stdin.") && !in1.equals("stdin")){
				File f=new File(in1);
				if(!f.exists() || !f.isFile()){
					in1=null;
//					throw new RuntimeException(in1+" does not exist.");
				}
			}
			if(in2!=null && !in2.contains(",")){
				File f=new File(in2);
				if(!f.exists() || !f.isFile()){
					in2=null;
//					throw new RuntimeException(in2+" does not exist.");
				}else if(in1.equalsIgnoreCase(in2)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		ReadWrite.MAX_ZIP_THREADS=Shared.threads()-1;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		
		boolean mm0set=false;
		KmerCountAbstract.minProb=0f;
		
		Parser parser=new Parser();
		parser.trimq=trimq;
		parser.minAvgQuality=minAvgQuality;
		parser.minReadLength=minReadLength;
		parser.maxReadLength=maxReadLength;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);}

			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("minentropy") || a.equals("entropy")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					minEntropyScore=Integer.parseInt(b);
				}else{
					useEntropy=Tools.parseBoolean(b);
				}
			}else if(a.equals("minoverlappingbases") || a.equals("minoverlapbases") || a.equals("minoverlap")){
				MIN_OVERLAPPING_BASES=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases0") || a.equals("minoverlapbases0") || a.equals("minoverlap0")){
				MIN_OVERLAPPING_BASES_0=Integer.parseInt(b);
			}else if(a.equals("minoverlapinsert") || a.equals("minoi")){
				MIN_OVERLAP_INSERT=Integer.parseInt(b);
			}else if(a.equals("minqo") || a.equals("minq")){
				MIN_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("maxq")){
				Read.MAX_MERGE_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("qualiters")){
				QUAL_ITERS=Tools.max(1, Integer.parseInt(b));
			}else if(a.equals("maxbadbases") || a.equals("maxbad") || a.equals("mismatches")){
				MAX_MISMATCHES=Integer.parseInt(b);
			}else if(a.equals("maxbadbases0") || a.equals("maxbad0") || a.equals("mismatches0")){
				MAX_MISMATCHES0=Integer.parseInt(b);
				mm0set=true;
			}else if(a.equals("margin")){
				MISMATCH_MARGIN=Integer.parseInt(b);
			}else if(a.equals("usemapping")){
				USE_MAPPING=Tools.parseBoolean(b);
			}else if(a.equals("bin")){
				bin=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Shared.setThreads(Integer.parseInt(b));
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.equals("outgood") || a.equals("outmerged") || a.equals("outm") || a.equals("out")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood1") || a.equals("outmerged1") || a.equals("outm1") || a.equals("out1")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood2") || a.equals("outmerged2") || a.equals("outm2") || a.equals("out2")){
				out2=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb") || a.equals("outu") || a.equals("outunmerged") || a.equals("outbad")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb1") || a.equals("outu1") || a.equals("outunmerged1") || a.equals("outbad1")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb2") || a.equals("outu2") || a.equals("outunmerged2") || a.equals("outbad2")){
				outb2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outinsert") || a.startsWith("outlength")){
				outinsert=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist") || a.equals("hist") || a.equals("histogram") || a.equals("ihist")){
				outhist=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist2") || a.equals("hist2")){
				outhist2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist3") || a.equals("hist3")){
				outhist3=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outputfailed")){
				OUTPUT_FAILED=Tools.parseBoolean(b);
			}else if(a.equals("mix")){
				MIX_BAD_AND_GOOD=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				NONZERO_ONLY=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("join") || a.equals("merge")){
				join=Tools.parseBoolean(b);
			}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
				trimByOverlap=Tools.parseBoolean(b);
			}else if(a.equals("useoverlap") || a.equals("usebases") || a.equals("matebyoverlap") || a.equals("matebybases")){
				MATE_BY_OVERLAP=Tools.parseBoolean(b);
			}else if(a.startsWith("skipmated")){
				SKIP_MATED_READS=Tools.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("trimonfailure") || a.equals("tof")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					TRIM_ON_OVERLAP_FAILURE=Integer.parseInt(b);
				}else{
					TRIM_ON_OVERLAP_FAILURE=(Tools.parseBoolean(b) ? 1 : 0);
				}
			}else if(a.equals("maxExpectedErrors") || a.equals("mee")){
				maxExpectedErrors=Float.parseFloat(b);
				//Read.AVERAGE_QUALITY_BY_PROBABILITY=false;
			}else if(a.equals("mi") || a.equals("minins") || a.equals("mininsert")){
				minInsert=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				KmerCountAbstract.minProb=Float.parseFloat(b);
				assert(KmerCountAbstract.minProb<1) : "minprob must be less than 1.  At 1, even kmers with 100% probablity of being error-free will be discarded.";
			}else if(a.equals("minqf") || a.equals("minqfilter")){
				filterMinq=(byte)Integer.parseInt(b);
			}else if(a.equals("k")){
				filterK=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				filterHashes=Integer.parseInt(b);
			}else if(a.equals("efilter")){
				useEfilter=Tools.parseBoolean(b);
			}else if(a.equals("kfilter")){
				useKFilter=Tools.parseBoolean(b);
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
			qtrim=((qtrimLeft||qtrimRight)&&trimq>=0);
//			minLenFraction=parser.minLenFraction;
			minAvgQuality=parser.minAvgQuality;
			minReadLength=parser.minReadLength;
			maxReadLength=(parser.maxReadLength<0 ? Integer.MAX_VALUE : parser.maxReadLength);
//			maxNs=parser.maxNs;
			untrim=parser.untrim;
//			minGC=parser.minGC;
//			maxGC=parser.maxGC;
//			filterGC=parser.filterGC;
//			minTrimLength=(parser.minTrimLength>=0 ? parser.minTrimLength : minTrimLength);
//			requireBothBad=parser.requireBothBad;
		}
		
		BBMergeOverlapper.USE_MAPPING=USE_MAPPING;
		BBMergeOverlapper.MIN_OVERLAP_INSERT=Tools.min(MIN_OVERLAP_INSERT, minInsert);
		if(verbose){BBMergeOverlapper.verbose=true;}
		
		if(trimByOverlap){
			join=false;
		}
		
		if(!mm0set){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
		}
		
		if(MAX_MISMATCHES0<MAX_MISMATCHES){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
			System.err.println("MAX_MISMATCHES0 was set to "+MAX_MISMATCHES0+" to remain >=MAX_MISMATCHES");
		}
		
		if(MISMATCH_MARGIN>MAX_MISMATCHES){
			MISMATCH_MARGIN=MAX_MISMATCHES;
			System.err.println("MISMATCH_MARGIN was set to "+MISMATCH_MARGIN+" to remain >=MAX_MISMATCHES");
		}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		if(in2==null && in1!=null && in1.contains("#") && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		
		if(out2==null && out1!=null && out1.contains("#")){
			out2=out1.replaceFirst("#", "2");
			out1=out1.replaceFirst("#", "1");
		}
		
		if(outb2==null && outb1!=null && outb1.contains("#")){
			outb2=outb1.replaceFirst("#", "2");
			outb1=outb1.replaceFirst("#", "1");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outb1, outb2, outinsert, outhist, outhist2, outhist3)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out1+", "+out2+", "+outb1+", "+outb2+", "+outinsert+", "+outhist+", "+outhist2+", "+outhist3+"\n");
		}
		
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
		
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		if(THREADS<1){THREADS=Shared.threads();}
		
		if(useKFilter){
			
			long memory=Runtime.getRuntime().maxMemory();
			long tmemory=Runtime.getRuntime().totalMemory();
			
			final long usable=(long)Tools.max(((memory-100000000-THREADS*8000000)*.70), memory*0.40);
			long mem=usable;
			
//			assert(false) : memory+", "+tmemory+", "+usable;
			
			filterCells=(mem*8)/filterBits;
			if(filterCells<80000000){
				System.err.println("Disabling filter due to low memory.  To enable the filter, increase the -Xmx flag.");
				useKFilter=false;
			}
		}
		useMEEfilter=maxExpectedErrors>0;
	}
	
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
		
		if(useKFilter){
			KmerCountAbstract.CANONICAL=true;
//			assert(false) : filterCells;
			filter=KmerCount7MTA.makeKca(in1, in2, extraFiles, filterK, filterBits, 0, filterCells, filterHashes, filterMinq, true, filterReads, 1, 1, 1, 1, null);
			System.err.println("Made prefilter:   \t"+filter.toShortString(filterHashes));
			double uf=filter.usedFraction();
			if(uf>0.7){
				System.err.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.98 ? "incredibly" : uf>0.95 ? "extremely" : uf>0.9 ? "very" : 
					uf>0.8 ? "fairly" : "somewhat")+" full, which may increase false-positive kmers.  Ideal load is under 70% used." +
					"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; " +
						"or increase the values of the minprob flag to reduce spurious kmers.");
			}
		}
		
		runPhase(join, maxReads, false);
		
		double stdev=0;
		if(histTotal!=null){
			stdev=Tools.standardDeviationHistogram(histTotal);
		}
		
		if(outhist!=null){
			StringBuilder sb=new StringBuilder();
			
			sb.append("#Mean\t"+String.format("%.3f", Tools.averageHistogram(histTotal))+"\n");
			sb.append("#Median\t"+Tools.percentile(histTotal, 0.5)+"\n");
			sb.append("#Mode\t"+Tools.calcMode(histTotal)+"\n");
			sb.append("#STDev\t"+String.format("%.3f", Tools.standardDeviationHistogram(histTotal))+"\n");
			sb.append("#InsertSize\tCount\n");
			for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				if(x>0 || !NONZERO_ONLY){
					sb.append(i+"\t"+x+"\n");
				}
			}
			ReadWrite.writeStringInThread(sb, outhist);
		}
		
		if(outhist2!=null){
			StringBuilder sb=new StringBuilder();
			sb.append("#InsertSize\tCount\n");
			int start=(histTotal.length==0 || histTotal[0]>0 ? 0 : 1);
			for(int i=start; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(x+"\n");
			}
			ReadWrite.writeStringInThread(sb, outhist2);
		}

		if(outhist3!=null){
			
			if(!new File(outhist3).exists()){
				StringBuilder sb=new StringBuilder();
				for(int i=0; i<histTotal.length; i+=bin){
					sb.append(i+"\n");
				}
				ReadWrite.writeString(sb, outhist3);
			}
			
			StringBuilder sb=new StringBuilder();
			TextFile tf=new TextFile(outhist3, false, false);
			for(int i=0; i<histTotal.length; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length && i<=insertMaxTotal; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(tf.readLine()+"\t"+x+"\n");
			}
			tf.close();
			ReadWrite.writeStringInThread(sb, outhist3);
		}
		
		
		ttotal.stop();
		System.err.println("Total time: "+ttotal+"\n");
		
		long sum=correctCountTotal+incorrectCountTotal;
		
		double div=100d/readsProcessedTotal;
		System.err.println("Pairs:       \t"+readsProcessedTotal);
		System.err.println("Joined:      \t"+sum+String.format((sum<10000 ? "       " : "   ")+"\t%.3f%%", sum*div));
		if(FASTQ.PARSE_CUSTOM){
			System.err.println("Correct:     \t"+correctCountTotal+String.format((correctCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", correctCountTotal*div));
			System.err.println("Incorrect:   \t"+incorrectCountTotal+String.format((incorrectCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", incorrectCountTotal*div));
		}
		System.err.println("Ambiguous:   \t"+ambiguousCountTotal+String.format((ambiguousCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", ambiguousCountTotal*div));
		System.err.println("No Solution: \t"+noSolutionCountTotal+String.format((noSolutionCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", noSolutionCountTotal*div));
		if(minInsert>0){System.err.println("Too Short:   \t"+tooShortCountTotal+String.format((tooShortCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooShortCountTotal*div));}
		if(maxReadLength<Integer.MAX_VALUE){System.err.println("Too Long:    \t"+tooLongCountTotal+String.format((tooLongCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooLongCountTotal*div));}
		System.err.println("Avg Insert:          \t\t"+String.format("%.1f", (insertSumCorrectTotal+insertSumIncorrectTotal)*1d/(correctCountTotal+incorrectCountTotal)));
		System.err.println("Standard Deviation:  \t\t"+String.format("%.1f", stdev));
		System.err.println("Mode:                \t"+Tools.calcMode(histTotal));
		if(FASTQ.PARSE_CUSTOM){
			System.err.println("Avg Insert Correct:  \t\t"+String.format("%.1f", (insertSumCorrectTotal)*1d/(correctCountTotal)));
			System.err.println("Avg Insert Incorrect:\t\t"+String.format("%.1f", (insertSumIncorrectTotal)*1d/(incorrectCountTotal)));
		}
		
		System.err.println();
		System.err.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
		System.err.println("90th percentile:     \t"+Tools.percentile(histTotal, .9));
		System.err.println("75th percentile:     \t"+Tools.percentile(histTotal, .75));
		System.err.println("50th percentile:     \t"+Tools.percentile(histTotal, .5));
		System.err.println("25th percentile:     \t"+Tools.percentile(histTotal, .25));
		System.err.println("10th percentile:     \t"+Tools.percentile(histTotal, .1));
	}
	
	public void runPhase(boolean join, long maxReads, boolean perfectonly){
		
		Timer talign=new Timer();
		
		ConcurrentReadOutputStream rosgood=null;
		ConcurrentReadOutputStream rosbad=null;
		ConcurrentReadOutputStream rosinsert=null;
		
		if(out1!=null){
			if(join==true){
				if(out2==null){System.err.println("Writing mergable reads merged.");}
				else{
					System.err.println("WARNING: 2 output files specified even though 'merge=true'.  out2 will be ignored.");
					out2=null;
				}
			}else{
				if(out2==null){System.err.println("Writing mergable reads interleaved.");}
				else{System.err.println("Writing mergable reads unmerged in two files.");}
			}

			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(8, 2*THREADS);
			rosgood=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosgood.start();
		}
		
		if(outb1!=null){

			final FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, append, true);
			final FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, append, true);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosbad=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final int buff=Tools.max(16, 2*THREADS);
			
			String out1=outinsert.replaceFirst("#", "1");

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			
			ReadStreamWriter.HEADER=header();
			final FileFormat ff=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, true);
			rosinsert=ConcurrentReadOutputStream.getStream(ff, null, null, null, buff, null, false);
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			System.err.println("Started output threads.");
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		assert(paired);//Fails on empty files.
		if(verbose){System.err.println("Paired: "+paired);}
		
		talign.start();
		
		
		MateThread[] pta=new MateThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new MateThread(cris, rosgood, rosbad, rosinsert, join, perfectonly, trimByOverlap);
			pta[i].start();
		}

		insertMinTotal=999999999;
		insertMaxTotal=0;
		
		readsProcessedTotal=0;
		matedCountTotal=0;
		correctCountTotal=0;
		ambiguousCountTotal=0;
		tooShortCountTotal=0;
		tooLongCountTotal=0;
		incorrectCountTotal=0;
		noSolutionCountTotal=0;
		insertSumCorrectTotal=0;
		insertSumIncorrectTotal=0;
		
		Arrays.fill(histTotal, 0);
		
		for(int i=0; i<pta.length; i++){
			MateThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				readsProcessedTotal+=ct.pairsProcessed;
				matedCountTotal+=ct.matedCount;
				correctCountTotal+=ct.correctCount;
				ambiguousCountTotal+=ct.ambiguousCount;
				tooShortCountTotal+=ct.tooShortCount;
				tooLongCountTotal+=ct.tooLongCount;
				incorrectCountTotal+=ct.incorrectCount;
				noSolutionCountTotal+=ct.noSolutionCount;
				insertSumCorrectTotal+=ct.insertSumCorrect;
				insertSumIncorrectTotal+=ct.insertSumIncorrect;
				
				basesTrimmedTotal+=ct.basesTrimmedT;
				readsTrimmedTotal+=ct.readsTrimmedT;

				insertMinTotal=Tools.min(ct.insertMin, insertMinTotal);
				insertMaxTotal=Tools.max(ct.insertMax, insertMaxTotal);
				
//				System.err.println(ct.insertMin+", "+ct.insertMax);
				
				if(ct.hist!=null){
					for(int h=0; h<ct.hist.length; h++){
						histTotal[h]+=ct.hist[h];
					}
				}
			}
		}
		
		System.err.println("Finished reading");
		errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad, rosinsert);
		
		talign.stop();
//		System.err.println("Align time: "+talign);
	}
	


	public static String header(){
		return "#id\tnumericID\tinsert\tstatus\thashHits\thashMisses\tscore\tsum\tvotes\n";
	}
	
	
	private class MateThread extends Thread{
		
		
		public MateThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosgood_, ConcurrentReadOutputStream rosbad_, ConcurrentReadOutputStream rosinsert_,
				boolean joinReads_, boolean joinperfectonly_, boolean trimByOverlap_) {
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosinsert=rosinsert_;
			joinReads=joinReads_;
			trimReadsByOverlap=trimByOverlap_;
			joinperfectonly=joinperfectonly_;
			
			if(useEntropy){
				kmerCounts=new short[1<<(2*entropyK)];
			}else{
				kmerCounts=null;
			}
		}
		
		
		@Override
		public void run(){
			processMate();
		}

		private void processMate() {
			
			final boolean USE_MAPPING=BBMergeOverlapper.USE_MAPPING;
			final boolean ignoreMappingStrand=BBMergeOverlapper.IGNORE_MAPPING_STRAND;
			final int[] rvector=new int[5];
			assert(USE_MAPPING || MATE_BY_OVERLAP);

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null);
			}


			while(reads!=null && reads.size()>0){

				ArrayList<Read> listg=(rosgood==null ? null : new ArrayList<Read>());
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>());
				ArrayList<Read> listi=(rosinsert==null ? null : new ArrayList<Read>());

				for(Read r1 : reads){
					final Read r2=r1.mate;
					
					TrimRead tr1=null, tr2=null;
					
					boolean remove=false;
					if(qtrim){
						if(untrim){
							if(r1!=null){
								tr1=TrimRead.trim(r1, qtrimLeft, qtrimRight, trimq, 1);
								int x=(tr1==null ? 0 : tr1.leftTrimmed+tr1.rightTrimmed);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2!=null){
								tr2=TrimRead.trim(r2, qtrimLeft, qtrimRight, trimq, 1);
								int x=(tr2==null ? 0 : tr2.leftTrimmed+tr2.rightTrimmed);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
						}else{
							if(r1!=null){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2!=null){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
						}
					}

					if(!remove){
						int rlen=(r1==null ? 0 : r1.length());
						int rlen2=(r1.mateLength());
						if(minReadLength>0){
							remove&=(rlen<minReadLength && rlen2<minReadLength);
						}
						if(r1.quality!=null || r2.quality!=null){
							if(minAvgQuality>0){
								remove=remove && (r1.avgQuality()<minAvgQuality || r2.avgQuality()<minAvgQuality);
							}
							if(useMEEfilter){
								remove=remove && (r1.expectedErrors()>maxExpectedErrors || r2.expectedErrors()>maxExpectedErrors);
							}
						}
						if(remove){
							basesTrimmedT+=(rlen+rlen2);
							remove=true;
						}
					}

					if(!remove){
						int trueSize=-1;
						
						if(r2!=null){
							if(FASTQ.PARSE_CUSTOM){
								r1.setMapped(true);
								r2.setMapped(true);
								trueSize=Read.insertSizeMapped(r1, r2, ignoreMappingStrand);
								if(verbose){System.err.println("True Insert: "+trueSize);}
							}
							r2.reverseComplement();
						}
						
						pairsProcessed++;
						
//						Arrays.fill(rvector, 0);

						int bestInsert=-1;
						int bestScore=-1;
						int bestGood=-1;
						int bestBad=999999;
						boolean ambig, tooShort=false, tooLong=false;

						//					assert(false) : r+"\n"+(USE_MAPPING)+", "+(r.chrom==r.mate.chrom)+", "+()+", "+()+", "+()+", "+()+", ";

						if(r2==null){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=r1.length();
							assert(r1.length()==r1.insert()) : r1.length()+" != "+r1.insert()+"; actual = "+trueSize;
							//						if(bestInsert!=trueSize){
							//							System.err.println("Bad insert size for pre-joined read "+r.numericID+": len="+r.length()+", insert="+r.insert()+", actual="+trueSize);
							//						}
							ambig=false;
						}else if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && ((r1.mapped() || r1.synthetic()) && (r2.mapped() || r2.synthetic()))){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=trueSize;
							ambig=false;
						}else if(SKIP_MATED_READS && r1.insertvalid() && r1.insert()>0){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=r1.insert();
							ambig=false;
						}else{
							if(MATE_BY_OVERLAP){
								final int len1=r1.length(), len2=r2.length();
								assert(QUAL_ITERS>0);
								ambig=false;
								
								int minOverlap=MIN_OVERLAPPING_BASES;
								if(useEntropy){
									int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, entropyK, kmerCounts, minEntropyScore);
									int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, entropyK, kmerCounts, minEntropyScore);
									minOverlap=Tools.max(minOverlap, Tools.max(a, b));
								}
								
								for(int i=0; i<QUAL_ITERS && bestInsert<0 /*&& !ambig*/; i++){
									
									int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-i, minOverlap+i, 
											MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, (byte)(MIN_QUALITY-2*i));
									if(x>-1){
										//											System.err.println(trims);
										bestInsert=x;
										bestScore=rvector[0];
										bestGood=rvector[1];
										bestBad=rvector[2];
										ambig=(rvector[4]==1);
//										if(!ambig){break;}
										break;
									}
									
//									if(x<=0 || ambig){
//										x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-i+20, minOverlap+i+20, 
//												MISMATCH_MARGIN, MAX_MISMATCHES0+1, MAX_MISMATCHES+1, (byte)(MIN_QUALITY-2*i));
//										if(x>-1 && (rvector[4]!=1)){
//											//											System.err.println(trims);
//											bestInsert=x;
//											bestScore=rvector[0];
//											bestGood=rvector[1];
//											bestBad=rvector[2];
//											ambig=(rvector[4]==1);
//											break;
//										}
//									}
									
								}
								
								
								if(loose && bestInsert<0){//TODO check for estimated number of overlap errors
									int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap+2, 
											MISMATCH_MARGIN, MAX_MISMATCHES0+1, MAX_MISMATCHES+1, MIN_QUALITY-1);
									if(x>-1){
										bestInsert=x;
										bestScore=rvector[0];
										bestGood=rvector[1];
										bestBad=rvector[2];
										ambig=(rvector[4]==1);
									}
								}
//								if(bestInsert<0 && MISMATCH_MARGIN>1){
//									int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-1, MIN_OVERLAPPING_BASES+2, 
//											MISMATCH_MARGIN-1, MAX_MISMATCHES, MIN_QUALITY-1);
//									if(x>-1){
//										bestInsert=x;
//										bestScore=rvector[0];
//										bestGood=rvector[1];
//										bestBad=rvector[2];
//										ambig=(rvector[4]==1);
//									}
//								}
//								if(bestInsert<0){
//									int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-1, MIN_OVERLAPPING_BASES+2, 
//											MISMATCH_MARGIN-1, MAX_MISMATCHES+1, MIN_QUALITY);
//									if(x>-1){
//										bestInsert=x;
//										bestScore=rvector[0];
//										bestGood=rvector[1];
//										bestBad=rvector[2];
//										ambig=(rvector[4]==1);
//									}
//								}
								
								//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
								if(loose && useEfilter && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
									float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
									if((bestExpected+.05f)*10f<bestBad){ambig=true;}
								}
								
								for(int trims=0, q=trimq; trims<TRIM_ON_OVERLAP_FAILURE && !qtrim && bestInsert<0 /*&& !ambig*/; trims++, q+=8){
//									System.err.println(trims+", "+q);
									Serializable old1=r1.obj;
									Serializable old2=r2.obj;
									tr1=TrimRead.trim(r1, false, true, q, 1+len1*4/10); //r1.length());
									tr2=TrimRead.trim(r2, true, false, q, 1+len2*4/10); //r2.length());
									r1.obj=old1;
									r2.obj=old2;
									if(tr1!=null || tr2!=null){
//										System.err.println(r1.length()+", "+r2.length());
										int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap, MISMATCH_MARGIN, 
												MAX_MISMATCHES0, MAX_MISMATCHES, MIN_QUALITY);
										if(x>-1){
//											System.err.println(trims);
											bestInsert=x;
											bestScore=rvector[0];
											bestGood=rvector[1];
											bestBad=rvector[2];
											ambig=(rvector[4]==1);
											trims=TRIM_ON_OVERLAP_FAILURE;
										}else{
											if(tr1!=null){tr1.untrim();}
											if(tr2!=null){tr2.untrim();}
										}
									}
								}
								
//								if(!ambig && bestInsert<0 && len1>=100 && len2>=100){
//									int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0+20, minOverlap+20, 
//											MISMATCH_MARGIN, MAX_MISMATCHES0+1, MAX_MISMATCHES+1, (byte)(MIN_QUALITY));
//									if(x>-1){
//										//											System.err.println(trims);
//										bestInsert=x;
//										bestScore=rvector[0];
//										bestGood=rvector[1];
//										bestBad=rvector[2];
//										ambig=(rvector[4]==1);
//									}
//								}
								
							}else{
								ambig=false;
								bestInsert=-1;
							}
						}

						tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
						tooLong=(!ambig && !tooShort && bestInsert>0 && bestInsert>maxReadLength);
						
						if(joinperfectonly && bestBad>0){ambig=true;}
						
						Read joined=null;
						if(bestInsert>-1 && !ambig && !tooShort && r2!=null && (joinReads || useKFilter)){
							joined=r1.joinRead(bestInsert);
							if(useKFilter){
								int cov=BBMergeOverlapper.minCoverage(joined, filter, filterK, true, filterCutoff);
								if(cov<filterCutoff){bestInsert=-1;}
							}
							if(!joinReads || bestInsert<0){joined=null;}
						}

						if(ambig){ambiguousCount++;}
						else if(tooShort){tooShortCount++;}
						else if(tooLong){tooLongCount++;}
						else if(bestInsert==-1){noSolutionCount++;}
						else if(bestInsert==trueSize){correctCount++;insertSumCorrect+=bestInsert;}
						else{incorrectCount++;insertSumIncorrect+=bestInsert;}

						if(bestInsert>-1 && !ambig){
							insertMin=Tools.min(bestInsert, insertMin);
							insertMax=Tools.max(bestInsert, insertMax);
							hist[Tools.min(bestInsert, hist.length-1)]++;
						}
						r1.setInsert(ambig ? -1 : bestInsert);
						
						if(OUTPUT_FAILED || bestInsert>-1){

							if(untrim && (ambig || bestInsert<0 || !joinReads)){
								if(tr1!=null){tr1.untrim();}
								if(tr2!=null){tr2.untrim();}
							}
							
							if((ambig || bestInsert<0 || tooShort) && (rosbad!=null || !MIX_BAD_AND_GOOD)){
								if(rosbad!=null){
									if(listb!=null){listb.add(r1);}
								}
							}else{
								if(listg!=null){
									if(joined==null){
										joined=r1;
										if(r2!=null){
											if(joinReads){joined=r1.joinRead(bestInsert);}
											else if(trimReadsByOverlap){
												int trimLim=bestInsert-1;
												if(trimLim<r1.length()){
													if(verbose){System.err.println("Overlap right trimming r1 to "+0+", "+(trimLim));}
													int x=TrimRead.trimToPosition(r1, 0, trimLim, 1);
													if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r1.bases));}
												}
												if(trimLim<r2.length()){
													if(verbose){System.err.println("Overlap right trimming r2 to "+0+", "+(trimLim));}
													int x=TrimRead.trimToPosition(r2, 0, trimLim, 1);
													if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r2.bases));}
												}
											}
										}
									}
									listg.add(joined);
								}
							}

							if(rosinsert!=null){
								StringBuilder sb=new StringBuilder(40);
								sb.append(r1.id==null ? r1.numericID+"" : r1.id).append('\t');
								sb.append(r1.numericID).append('\t');

								sb.append(bestInsert);
								sb.append('\t');

								if(bestInsert<0){sb.append('F');}//Failed
								else if(ambig){sb.append('A');} //Ambiguous
								else if(tooShort){sb.append('S');} //Short
								else if(bestInsert>0 && bestBad<1){sb.append('P');} //Perfect
								else{sb.append('I');}//Imperfect

								if(bestInsert>0){
									sb.append("\t"+bestGood+"\t"+bestBad+"\t"+bestScore);
								}
								r1.obj=sb;
								listi.add(r1);
							}
						}
						if(r2!=null){r2.reverseComplement();}
					}
				}

				if(rosgood!=null){rosgood.add(listg, ln.id);}
				if(rosbad!=null){rosbad.add(listb, ln.id);}
				if(rosinsert!=null){rosinsert.add(listi, ln.id);}

				//			System.err.println("returning list");
				cris.returnList(ln.id, ln.list.isEmpty());
				//			System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			System.err.println("reads: "+(reads==null ? "null" : reads.size()));
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}
		
		/*--------------------------------------------------------------*/

		final long[] hist=new long[histlen];
		final short[] kmerCounts;

		long pairsProcessed=0;
		long matedCount=0;
		long correctCount=0;
		long ambiguousCount=0;
		long tooShortCount=0;
		long tooLongCount=0;
		long incorrectCount=0;
		long noSolutionCount=0;
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		int insertMax=0;
		int insertMin=999999999;
		
		long basesTrimmedT=0;
		long readsTrimmedT=0;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream rosgood;
		private final ConcurrentReadOutputStream rosbad;
		private final ConcurrentReadOutputStream rosinsert;
		
		private final boolean joinReads;
		private final boolean joinperfectonly;
		private final boolean trimReadsByOverlap;
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1;
	private String in2;

	private String out1=null;
	private String out2=null;
	private String outb1=null;
	private String outb2=null;
	private String outinsert=null;
	private String outhist=null;
	private String outhist2=null;
	private String outhist3=null;
	
	private long maxReads=-1;
	private boolean join=true;
	private boolean trimByOverlap=false;

	private boolean useEfilter=false;
	private boolean useMEEfilter=false;
	
	private boolean useKFilter=false;
	private int filterBits=1;
	private int filterCutoff=(1<<filterBits)-1;
	private long filterCells=-1;
	private int filterK=31;
	private int filterHashes=1;
	private long filterReads=-1;
	private KCountArray filter;
	private int filterMinq=2;
	private ArrayList<String> extraFiles;

	private boolean useEntropy=true;
	private int entropyK=3;
	private int minEntropyScore=39;//30 loose;//39 normal;//44 strict;
	
	private byte qin=-1;
	private byte qout=-1;
	
	/*--------------------------------------------------------------*/
	
	static boolean errorState=false;
	
	static boolean qtrimRight=false;
	static boolean qtrimLeft=false;
	static boolean untrim=false;
	static byte trimq=6;
	static byte minAvgQuality=0;
	static float maxExpectedErrors=0;
	static int minReadLength=0;
	static int maxReadLength=-1;
	static int minInsert=35;
	static boolean qtrim=false;
	static int TRIM_ON_OVERLAP_FAILURE=1;
	static int QUAL_ITERS=3;
	
	static boolean strict=false;
	static boolean vstrict=false;
	static boolean loose=false;
	static boolean vloose=false;
	static boolean fast=false;
	
	private static final int histlen=2000; 
	static long[] histTotal=new long[histlen];
	static int bin=1;

	static long readsProcessedTotal=0;
	static long matedCountTotal=0;
	static long correctCountTotal=0;
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long tooLongCountTotal=0;
	static long incorrectCountTotal=0;
	static long noSolutionCountTotal=0;
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	static long basesTrimmedTotal=0;
	static long readsTrimmedTotal=0;
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;
	
	private static int MIN_OVERLAPPING_BASES=11;
	private static int MIN_OVERLAPPING_BASES_0=8;
	private static int MISMATCH_MARGIN=2;

	public static int MAX_MISMATCHES=3;
	public static int MAX_MISMATCHES0=3;
	public static byte MIN_QUALITY=10;
	
	public static boolean USE_MAPPING=false;
	public static int MIN_OVERLAP_INSERT=25;
	
	private static boolean MATE_BY_OVERLAP=true;
	private static boolean SKIP_MATED_READS=false;
	private static boolean OUTPUT_FAILED=true;
	private static boolean MIX_BAD_AND_GOOD=false;
	private static boolean NONZERO_ONLY=true;
	
	private static boolean overwrite=true;
	private static boolean append=false;
	private static boolean verbose=false;
	
	private static int THREADS=-1;
	private static float version=5.5f;
	
}
