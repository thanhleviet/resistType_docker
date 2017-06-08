package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;

import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;

import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;

/**
 * @author Brian Bushnell
 * @date Sep 11, 2012
 *
 */
public class ReformatReads {

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		ReformatReads rr=new ReformatReads(args);
		rr.process(t);
	}
	
	public ReformatReads(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		SamLine.SET_FROM_OK=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens
			
			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("null") || a.equals(parser.in2)){
				// do nothing
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("sample") || a.equals("samplereads") || a.equals("samplereadstarget") || a.equals("srt")){
				sampleReadsTarget=Tools.parseKMG(b);
				sampleReadsExact=(sampleReadsTarget>0);
			}else if(a.equals("samplebases") || a.equals("samplebasestarget") || a.equals("sbt")){
				sampleBasesTarget=Tools.parseKMG(b);
				sampleBasesExact=(sampleBasesTarget>0);
			}else if(a.equals("addslash")){
				addslash=Tools.parseBoolean(b);
			}else if(a.equals("uniquenames")){
				uniqueNames=Tools.parseBoolean(b);
			}else if(a.equals("verifyinterleaved") || a.equals("verifyinterleaving") || a.equals("vint")){
				verifyinterleaving=Tools.parseBoolean(b);
			}else if(a.equals("verifypaired") || a.equals("verifypairing") || a.equals("vpair")){
				verifypairing=Tools.parseBoolean(b);
			}else if(a.equals("allowidenticalnames") || a.equals("ain")){
				allowIdenticalPairNames=Tools.parseBoolean(b);
			}else if(a.equals("rcompmate") || a.equals("rcm")){
				reverseComplimentMate=Tools.parseBoolean(b);
				outstream.println("Set RCOMPMATE to "+reverseComplimentMate);
			}else if(a.equals("rcomp") || a.equals("rc")){
				reverseCompliment=Tools.parseBoolean(b);
				outstream.println("Set RCOMP to "+reverseCompliment);
			}else if(a.equals("deleteempty") || a.equals("deletempty") || a.equals("delempty") || a.equals("def")){
				deleteEmptyFiles=Tools.parseBoolean(b);
			}else if(a.equals("mappedonly")){
				mappedOnly=Tools.parseBoolean(b);
			}else if(a.equals("unmappedonly")){
				unmappedOnly=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					parser.in1=b.replace("#", "1");
					parser.in2=b.replace("#", "2");
				}
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Download parser fields
			
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			testsize=parser.testsize;
			trimBadSequence=parser.trimBadSequence;
			breakLength=parser.breakLength;
			
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			minAvgQuality=parser.minAvgQuality;
			chastityFilter=parser.chastityFilter;
			maxNs=parser.maxNs;
			minReadLength=parser.minReadLength;
			maxReadLength=parser.maxReadLength;
			minLenFraction=parser.minLenFraction;
			requireBothBad=parser.requireBothBad;
			minGC=parser.minGC;
			maxGC=parser.maxGC;
			filterGC=parser.filterGC;

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			outsingle=parser.outsingle;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(ReadStats.INDEL_HIST_FILE!=null){SamLine.CONVERT_CIGAR_TO_MATCH=true;}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		if(SamLine.setxs && !SamLine.setintron){SamLine.INTRON_LIMIT=10;}
		qtrim=qtrimLeft||qtrimRight;
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(verifyinterleaving || (verifypairing && in2==null)){
			verifypairing=true;
			setInterleaved=true;
//			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
			if(!parser.setOut){
				System.err.println("No output stream specified.  To write to stdout, please specify 'out=stdout.fq' or similar.");
//				out1="stdout";
			}
		}
		
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		if(outsingle!=null && outsingle.equalsIgnoreCase("null")){outsingle=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outsingle)){
			System.err.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outsingle) || !ReadStats.testFiles(false)){
			throw new RuntimeException("Duplicate filenames are not allowed.");
		}
		
		FASTQ.PARSE_CUSTOM=parsecustom;
		
		{
			byte qin=Parser.qin, qout=Parser.qout;
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
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutsingle=FileFormat.testOutput(outsingle, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		assert(ReadStats.testFiles(true)) : "Existing output files specified, but overwrite==false";
		assert(ReadStats.testFiles(false)) : "Duplicate or output files specified";

//		System.err.println("\n"+ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1+"\n");
//		assert(false) : ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1;

		if(ffin1!=null && ffout1!=null && ffin1.samOrBam()){
			if(ffout1.samOrBam()){
				useSharedHeader=true;
				SamLine.CONVERT_CIGAR_TO_MATCH=true;
			}else if(ffout1.bread()){
				SamLine.CONVERT_CIGAR_TO_MATCH=true;
			}
		}

		nameMap1=(uniqueNames ? new HashMap<String, Integer>() : null);
		nameMap2=(uniqueNames ? new HashMap<String, Integer>() : null);
	}
	
	void process(Timer t){
		
		long readsRemaining=0;
		long basesRemaining=0;
		
		if(sampleReadsExact || sampleBasesExact){
			long[] counts=countReads(in1, in2, maxReads);
			readsRemaining=counts[0];
			basesRemaining=counts[2];
			setSampleSeed(sampleseed);
		}
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, colorspace, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		if(verbose){
			System.err.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
//		}
		
		assert(!paired || breakLength<1) : "Paired input cannot be broken with 'breaklength'";

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader);
			ros.start();
		}else{ros=null;}
		
		final ConcurrentReadOutputStream rosb;
		if(outsingle!=null){
			final int buff=4;
			
			rosb=ConcurrentReadOutputStream.getStream(ffoutsingle, null, buff, null, useSharedHeader);
			rosb.start();
		}else{rosb=null;}
		final boolean discardTogether=(!paired || (outsingle==null && !requireBothBad));
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		//Only used with deleteEmptyFiles flag
		long readsOut1=0;
		long readsOut2=0;
		long readsOutSingle=0;
		
		long basesOut1=0;
		long basesOut2=0;
		long basesOutSingle=0;
		
		long basesFTrimmedT=0;
		long readsFTrimmedT=0;
		
		long basesQTrimmedT=0;
		long readsQTrimmedT=0;
		
		long lowqBasesT=0;
		long lowqReadsT=0;
		
		long badGcBasesT=0;
		long badGcReadsT=0;
		
		long readShortDiscardsT=0;
		long baseShortDiscardsT=0;
		
		long unmappedReadsT=0;
		long unmappedBasesT=0;

		final boolean MAKE_QHIST=ReadStats.COLLECT_QUALITY_STATS;
		final boolean MAKE_QAHIST=ReadStats.COLLECT_QUALITY_ACCURACY;
		final boolean MAKE_MHIST=ReadStats.COLLECT_MATCH_STATS;
		final boolean MAKE_BHIST=ReadStats.COLLECT_BASE_STATS;
		
		final boolean MAKE_EHIST=ReadStats.COLLECT_ERROR_STATS;
		final boolean MAKE_INDELHIST=ReadStats.COLLECT_INDEL_STATS;
		final boolean MAKE_LHIST=ReadStats.COLLECT_LENGTH_STATS;
		final boolean MAKE_GCHIST=ReadStats.COLLECT_GC_STATS;
		final boolean MAKE_IDHIST=ReadStats.COLLECT_IDENTITY_STATS;
		
		final ReadStats readstats=(MAKE_QHIST || MAKE_MHIST || MAKE_BHIST || MAKE_QAHIST || MAKE_EHIST || MAKE_INDELHIST || MAKE_LHIST || MAKE_GCHIST || MAKE_IDHIST) ? 
				new ReadStats() : null;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(32));
				
				if(breakLength>0){
					breakReads(reads, breakLength, minReadLength);
				}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					if(readstats!=null){
						if(MAKE_QHIST){readstats.addToQualityHistogram(r1);}
						if(MAKE_BHIST){readstats.addToBaseHistogram(r1);}
						if(MAKE_MHIST){readstats.addToMatchHistogram(r1);}
						if(MAKE_QAHIST){readstats.addToQualityAccuracy(r1);}

						if(MAKE_EHIST){readstats.addToErrorHistogram(r1);}
						if(MAKE_INDELHIST){readstats.addToIndelHistogram(r1);}
						if(MAKE_LHIST){readstats.addToLengthHistogram(r1);}
						if(MAKE_GCHIST){readstats.addToGCHistogram(r1);}
						if(MAKE_IDHIST){readstats.addToIdentityHistogram(r1);}
					}
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
						if(reverseCompliment){r1.reverseComplement();}
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
						if(reverseCompliment || reverseComplimentMate){r2.reverseComplement();}
					}
					
					if(verifypairing){
						String s1=r1==null ? null : r1.id;
						String s2=r2==null ? null : r2.id;
						boolean b=FASTQ.testPairNames(s1, s2, allowIdenticalPairNames);
						if(!b){
							outstream.println("Names do not appear to be correctly paired.\n"+s1+"\n"+s2+"\n");
							ReadWrite.closeStreams(cris, ros);
							System.exit(1);
						}
					}
					
					if(trimBadSequence){//Experimental
						if(r1!=null){
							int x=TrimRead.trimBadSequence(r1);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null){
							int x=TrimRead.trimBadSequence(r2);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(chastityFilter){
						if(r1!=null && r1.failsChastity()){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && r2.failsChastity()){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(mappedOnly){
						if(r1!=null && !r1.discarded() && (!r1.mapped() || r1.bases==null || r1.secondary())){
							r1.setDiscarded(true);
							unmappedBasesT+=initialLength1;
							unmappedReadsT++;
						}
						if(r2!=null && !r2.discarded() && (!r2.mapped() || r2.bases==null || r2.secondary())){
							r2.setDiscarded(true);
							unmappedBasesT+=initialLength2;
							unmappedReadsT++;
						}
					}else if(unmappedOnly){
						if(r1!=null && !r1.discarded() && (!r1.mapped() || r1.bases==null || r1.secondary())){
							r1.setDiscarded(true);
							unmappedBasesT+=initialLength1;
							unmappedReadsT++;
						}
						if(r2!=null && !r2.discarded() && (!r2.mapped() || r2.bases==null || r2.secondary())){
							r2.setDiscarded(true);
							unmappedBasesT+=initialLength2;
							unmappedReadsT++;
						}
					}
					
					if(filterGC && (initialLength1>0 || initialLength2>0)){
						final float gc;
						if(r2==null){
							gc=r1.gc();
						}else{
							gc=(r1.gc()*initialLength1+r2.gc()*initialLength2)/(initialLength1+initialLength2);
						}
						if(gc<minGC || gc>maxGC){
							if(r1!=null && !r1.discarded()){
								r1.setDiscarded(true);
								badGcBasesT+=initialLength1;
								badGcReadsT++;
							}
							if(r2!=null && !r2.discarded()){
								r2.setDiscarded(true);
								badGcBasesT+=initialLength2;
								badGcReadsT++;
							}
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
					
					if(qtrim){
						if(r1!=null && !r1.discarded()){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null && !r2.discarded()){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(minAvgQuality>0){
						if(r1!=null && !r1.discarded() && r1.avgQuality()<minAvgQuality){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && !r2.discarded() && r2.avgQuality()<minAvgQuality){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(maxNs>=0){
						if(r1!=null && !r1.discarded() && r1.countUndefined()>maxNs){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && !r2.discarded() && r2.countUndefined()>maxNs){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(minlen1>0 || minlen2>0 || maxReadLength>0){
//						assert(false) : minlen1+", "+minlen2+", "+maxReadLength+", "+r1.length();
						if(r1!=null && !r1.discarded()){
							int rlen=r1.length();
							if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
								r1.setDiscarded(true);
								readShortDiscardsT++;
								baseShortDiscardsT+=rlen;
							}
						}
						if(r2!=null && !r2.discarded()){
							int rlen=r2.length();
							if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
								r2.setDiscarded(true);
								readShortDiscardsT++;
								baseShortDiscardsT+=rlen;
							}
						}
					}
					
					boolean remove=false;
					if(r2==null){
						remove=r1.discarded();
					}else{
						remove=requireBothBad ? (r1.discarded() && r2.discarded()) : (r1.discarded() || r2.discarded());
					}
					
					if(remove){reads.set(idx, null);}
					else{
						if(uniqueNames){
							{
								Integer v=nameMap1.get(r1.id);
								if(v==null){
									nameMap1.put(r1.id, 1);
								}else{
									v++;
									nameMap1.put(r1.id, v);
									r1.id=r1.id+"_"+v;
								}
							}
							if(r2!=null){
								Integer v=nameMap2.get(r2.id);
								if(v==null){
									nameMap2.put(r2.id, 1);
								}else{
									v++;
									nameMap2.put(r2.id, v);
									r2.id=r2.id+"_"+v;
								}
							}
						}
						if(addslash){
							if(r1.id==null){r1.id=" "+r1.numericID;}
							if(!r1.id.contains(" /1")){r1.id+=" /1";}
							if(r2!=null){
								if(r2.id==null){r2.id=" "+r2.numericID;}
								if(!r2.id.contains(" /2")){r2.id+=" /2";}
							}
						}
					}
					
					if(singles!=null){
						if(r1.discarded() || (r2!=null && r2.discarded())){
							if(!r1.discarded()){
								Read r=r1.clone();
								r.mate=null;
								r.setPairnum(0);
								singles.add(r);
							}else if(r2!=null && !r2.discarded()){
								Read r=r2.clone();
								r.mate=null;
								r.setPairnum(0);
								singles.add(r);
							}
						}
					}
				}
				
				final ArrayList<Read> listOut;
				
//				assert(false) : sampleReadsExact+", "+sampleBasesExact;
				if(sampleReadsExact || sampleBasesExact){
					listOut=new ArrayList<Read>();
					if(sampleReadsExact){
						for(Read r : reads){
							if(r!=null){
								assert(readsRemaining>0) : readsRemaining;
								double prob=sampleReadsTarget/(double)(readsRemaining);
//								System.err.println("sampleReadsTarget="+sampleReadsTarget+", readsRemaining="+readsRemaining+", prob="+prob);
								if(randy.nextDouble()<prob){
									listOut.add(r);
									sampleReadsTarget--;
								}
							}
							readsRemaining--;
						}
					}else if(sampleBasesExact){
						for(Read r : reads){
							if(r!=null){
								assert(basesRemaining>0) : basesRemaining;
								int bases=r.length()+(r.mate==null ? 0 : r.mateLength());
								double prob=sampleBasesTarget/(double)(basesRemaining);
								if(randy.nextDouble()<prob){
									listOut.add(r);
									sampleBasesTarget-=bases;
								}
								basesRemaining-=bases;
							}
						}
					}
				}else{
					listOut=reads;
				}
//				if(deleteEmptyFiles){
					for(Read r : listOut){
						if(r!=null){
							readsOut1++;
							basesOut1+=r.length();
							if(r.mate!=null){
								readsOut2++;
								basesOut2+=r.mateLength();
							}
						}
					}
					if(singles!=null){
						for(Read r : singles){
							if(r!=null){
								readsOutSingle++;
								basesOutSingle+=r.length();
							}
						}
					}
//				}
				if(ros!=null){ros.add(listOut, ln.id);}
				if(rosb!=null){rosb.add(singles, ln.id);}

				cris.returnList(ln.id, false);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
//				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				assert(ln.list.isEmpty());
				cris.returnList(ln.id, true);
			}
		}
		
		errorState|=ReadStats.writeAll(paired);
		
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
		
		if(deleteEmptyFiles){
			deleteEmpty(readsOut1, readsOut2, readsOutSingle);
		}
		
//		System.err.println(cris.errorState()+", "+(ros==null ? "null" : (ros.errorState()+", "+ros.finishedSuccessfully())));
//		if(ros!=null){
//			ReadStreamWriter rs1=ros.getRS1();
//			ReadStreamWriter rs2=ros.getRS2();
//			System.err.println(rs1==null ? "null" : rs1.finishedSuccessfully());
//			System.err.println(rs2==null ? "null" : rs2.finishedSuccessfully());
//		}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		final long rawReadsIn=cris.readsIn(), rawBasesIn=cris.basesIn();
		final double rmult=100.0/rawReadsIn, bmult=100.0/rawBasesIn;
		final double rpmult=100.0/readsProcessed, bpmult=100.0/basesProcessed;
		
		outstream.println("Input:                  \t"+cris.readsIn()+" reads          \t"+
				cris.basesIn()+" bases");
		if(samplerate!=1f){
			outstream.println("Processed:              \t"+readsProcessed+" reads          \t"+
					basesProcessed+" bases");
		}
		
		if(qtrim || trimBadSequence){
			outstream.println("QTrimmed:               \t"+readsQTrimmedT+" reads ("+String.format("%.2f",readsQTrimmedT*rpmult)+"%) \t"+
					basesQTrimmedT+" bases ("+String.format("%.2f",basesQTrimmedT*bpmult)+"%)");
		}
		if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimModulo>0){
			outstream.println("FTrimmed:               \t"+readsFTrimmedT+" reads ("+String.format("%.2f",readsFTrimmedT*rpmult)+"%) \t"+
					basesFTrimmedT+" bases ("+String.format("%.2f",basesFTrimmedT*bpmult)+"%)");
		}
		if(minReadLength>0 || maxReadLength>0){
			outstream.println("Short Read Discards:    \t"+readShortDiscardsT+" reads ("+String.format("%.2f",readShortDiscardsT*rpmult)+"%) \t"+
					baseShortDiscardsT+" bases ("+String.format("%.2f",baseShortDiscardsT*bpmult)+"%)");
		}
		if(minAvgQuality>0 || maxNs>=0 || chastityFilter){
			outstream.println("Low quality discards:   \t"+lowqReadsT+" reads ("+String.format("%.2f",lowqReadsT*rpmult)+"%) \t"+
					lowqBasesT+" bases ("+String.format("%.2f",lowqBasesT*bpmult)+"%)");
		}
		if(filterGC){
			outstream.println("GC content discards:    \t"+badGcReadsT+" reads ("+String.format("%.2f",badGcReadsT*rpmult)+"%) \t"+
					badGcBasesT+" bases ("+String.format("%.2f",badGcBasesT*bpmult)+"%)");
		}
		final long ro=readsOut1+readsOut2+readsOutSingle, bo=basesOut1+basesOut2+basesOutSingle;
		outstream.println("Output:                 \t"+ro+" reads ("+String.format("%.2f",ro*rmult)+"%) \t"+
				bo+" bases ("+String.format("%.2f",bo*bmult)+"%)");
		
		outstream.println("\nTime:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length())+
					(qfin1==null ? 0 : new File(qfin1).length())+(qfin2==null ? 0 : new File(qfin2).length()));//*passes
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+String.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(verifypairing){
			outstream.println("Names appear to be correctly paired.");
		}
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void deleteEmpty(long readsOut1, long readsOut2, long readsOutSingle){
		deleteEmpty(readsOut1, ffout1, qfout1);
		deleteEmpty(readsOut2, ffout2, qfout2);
		deleteEmpty(readsOutSingle, ffoutsingle, null);
	}
	
	private void deleteEmpty(long count, FileFormat ff, String qf){
		try {
			if(ff!=null && count<1){
				String s=ff.name();
				if(s!=null && !ff.stdio() && !ff.devnull()){
					File f=new File(ff.name());
					if(f.exists()){
						f.delete();
					}
				}
				if(qf!=null){
					File f=new File(qf);
					if(f.exists()){
						f.delete();
					}
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

	public static void breakReads(ArrayList<Read> list, final int max, int min){
		if(!containsReadsOutsideSizeRange(list, min, max)){return;}
		assert(max>0 || min>0) : "min or max read length must be positive.";
		assert(max<1 || max>=min) : "max read length must be at least min read length: "+max+"<"+min;
		min=Tools.max(0, min);
		
		ArrayList<Read> temp=new ArrayList<Read>(list.size()*2);
		for(Read r : list){
			if(r==null || r.bases==null){
				temp.add(r);
			}else if(r.length()<min){
				temp.add(null);
			}else if(max<1 || r.length()<=max){
				temp.add(r);
			}else{
				final byte[] bases=r.bases;
				final byte[] quals=r.quality;
				final String name=r.id;
				final int limit=bases.length-min;
				for(int num=1, start=0, stop=max; start<limit; num++, start+=max, stop+=max){
					if(verbose){
						System.err.println(bases.length+", "+start+", "+stop);
						if(quals!=null){System.err.println(quals.length+", "+start+", "+stop);}
					}
					stop=Tools.min(stop, bases.length);
					byte[] b2=Arrays.copyOfRange(bases, start, stop);
					byte[] q2=(quals==null ? null : Arrays.copyOfRange(quals, start, stop));
					String n2=name+"_"+num;
					Read r2=new Read(b2, -1, -1, -1, n2, q2, r.numericID, r.flags);
					r2.setMapped(false);
					temp.add(r2);
				}
			}
		}
		list.clear();
		list.ensureCapacity(temp.size());
//		list.addAll(temp);
		for(Read r : temp){
			if(r!=null){list.add(r);}
		}
	}
	
	private static boolean containsReadsAboveSize(ArrayList<Read> list, int size){
		for(Read r : list){
			if(r!=null && r.bases!=null){
				if(r.length()>size){
					assert(r.mate==null) : "Read of length "+r.length()+">"+size+". Paired input is incompatible with 'breaklength'";
					return true;
				}
			}
		}
		return false;
	}
	
	private static boolean containsReadsOutsideSizeRange(ArrayList<Read> list, int min, int max){
		for(Read r : list){
			if(r!=null && r.bases!=null){
				if((max>0 && r.length()>max) || r.length()<min){
					assert(r.mate==null) : "Read of length "+r.length()+" outside of range "+min+"-"+max+". Paired input is incompatible with 'breaklength'";
					return true;
				}
			}
		}
		return false;
	}
	
	
	private long[] countReads(String fname1, String fname2, long maxReads){
		{
			String x=fname1.toLowerCase();
			if((x.equals("stdin") || x.startsWith("stdin.")) && !new File(fname1).exists()){
				throw new RuntimeException("Can't precount reads from standard in, only from a file.");
			}
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, colorspace, false, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Counting Reads");}
			cris.start(); //4567
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		long count=0, count2=0, bases=0;
		
		while(reads!=null && reads.size()>0){
			count+=reads.size();
			for(Read r : reads){
				bases+=r.length();
				count2++;
				if(r.mate!=null){
					bases+=r.mateLength();
					count2++;
				}
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln.id, ln.list.isEmpty());
		errorState|=ReadWrite.closeStream(cris);
		return new long[] {count, count2, bases};
	}
	
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
		outstream.println("Other parameters and their defaults:\n");
		outstream.println("overwrite=false  \tOverwrites files that already exist");
		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
		outstream.println("fastawrap=80     \tLength of lines in fasta output");
		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	
	public void setSampleSeed(long seed){
		randy=new Random();
		if(seed>-1){randy.setSeed(seed);}
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;
	private String outsingle=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	/** Tracks names to ensure no duplicate names. */
	private final HashMap<String,Integer> nameMap1, nameMap2;
	private boolean uniqueNames=false;

	private boolean reverseComplimentMate=false;
	private boolean reverseCompliment=false;
	private boolean verifyinterleaving=false;
	private boolean verifypairing=false;
	private boolean allowIdenticalPairNames=true;
	private boolean trimBadSequence=false;
	private boolean chastityFilter=false;
	private boolean deleteEmptyFiles=false;
	private boolean mappedOnly=false;
	private boolean unmappedOnly=false;
	/** Add /1 and /2 to paired reads */
	private boolean addslash=false;

	private long maxReads=-1;
	private float samplerate=1f;
	private long sampleseed=-1;
	private boolean sampleReadsExact=false;
	private boolean sampleBasesExact=false;
	private long sampleReadsTarget=0;
	private long sampleBasesTarget=0;
	
	private boolean qtrimRight=false;
	private boolean qtrimLeft=false;
	private final int forceTrimLeft;
	private final int forceTrimRight;
	/** Trim right bases of the read modulo this value. 
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	private final int forceTrimModulo;
	private byte trimq=6;
	private byte minAvgQuality=0;
	private int maxNs=-1;
	private int breakLength=0;
	private int maxReadLength=0;
	private int minReadLength=0;
	private float minLenFraction=0;
	private float minGC=0;
	private float maxGC=1;
	private boolean filterGC=false;
	/** Toss pair only if both reads are shorter than limit */ 
	private boolean requireBothBad=false;
	
	private boolean useSharedHeader;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	private final FileFormat ffoutsingle;
	
	private final boolean qtrim;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	private boolean colorspace=false;
	private boolean parsecustom=false;
	private boolean testsize=false;
	
	private Random randy;
	
}
