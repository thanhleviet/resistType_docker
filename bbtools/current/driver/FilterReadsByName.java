package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

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
import fileIO.TextFile;

import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 8, 2014
 *
 */
public class FilterReadsByName {

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		FilterReadsByName mb=new FilterReadsByName(args);
		mb.process(t);
	}
	
	public FilterReadsByName(String[] args){
		
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
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("names")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.add(s);
					}
				}
			}else if(a.equals("include") || a.equals("retain")){
				exclude=!Tools.parseBoolean(b);
			}else if(a.equals("exclude") || a.equals("remove")){
				exclude=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					parser.in1=b.replace("#", "1");
					parser.in2=b.replace("#", "2");
				}
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{
			String[] x=names.toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				if(new File(s).exists()){
					
					int[] vector=FileFormat.testFormat(s, false, false);
					final int type=vector[0];
					
					TextFile tf=new TextFile(s);
					String[] lines=tf.toStringLines();
					tf.close();
					
					if(lines!=null && lines.length>0 && lines[0]!=null){
						String line0=lines[0];
						if(type==FileFormat.FASTQ && line0.startsWith("@")){
//							assert(false) : "A: "+line0;
							int mod=0;
							for(String s2 : lines){
								if((mod&3)==0){names.add(s2.substring(1));}
								mod++;
							}
						}else if(vector[0]==FileFormat.FASTA && line0.startsWith(">")){
//							assert(false) : "B: "+line0;
							for(String s2 : lines){
								if(s2.startsWith(">")){names.add(s2.substring(1));}
							}
//							for(String s2 : names){assert(false) : s2;}
						}else if(vector[0]==FileFormat.SAM && line0.contains("\t")){
//							assert(false) : "C: "+line0;
							for(String s2 : lines){
								names.add(new SamLine(s2).qname);
							}
						}else{
//							assert(false) : "D: "+line0;
							for(String s2 : lines){
								names.add(s2);
							}
						}
					}
				}else{
					names.add(s);
				}
			}
		}
		
		{//Download parser fields
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
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
//			if(!parser.setOut){
//				out1="stdout";
//			}
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
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
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

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, false, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		if(verbose){
			outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
//		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long readsOut=0;
		long basesOut=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			
			
			while(reads!=null && reads.size()>0){
				
				ArrayList<Read> retain=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					String name=r1.id, prefix=null;
					for(int x=0; x<name.length(); x++){
						char c=name.charAt(x);
						if((c=='/' && x<name.length()-1 && name.charAt(x+1)=='1') || Character.isWhitespace(c)){
							prefix=name.substring(0, x);
							break;
						}
					}
					
					boolean match=(names.contains(name) || prefix!=null && names.contains(prefix));
					
//					assert(false) : names.contains(name)+", "+name+", "+prefix+", "+exclude;
					
					if(match!=exclude){
						retain.add(r1);
						{
							readsOut++;
							basesOut+=initialLength1;
						}
						if(r2!=null){
							readsOut++;
							basesOut+=initialLength2;
						}
					}
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
				}
				
				final ArrayList<Read> listOut=retain;
				
				if(ros!=null){ros.add(listOut, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadStats.writeAll(paired);
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

//		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
//		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
//
//		while(rpstring.length()<8){rpstring=" "+rpstring;}
//		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:               "+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Processed:    "+readsOut);
		outstream.println("Bases Processed:    "+basesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
//		outstream.println("Syntax:\n");
//		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
//		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
//		outstream.println("Other parameters and their defaults:\n");
//		outstream.println("overwrite=false  \tOverwrites files that already exist");
//		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
//		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
//		outstream.println("fastawrap=80     \tLength of lines in fasta output");
//		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
//		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
//		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private boolean exclude=true;
	
	private HashSet<String> names=new HashSet<String>();
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
