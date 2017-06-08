package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.FileFormat;

import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Jan 12, 2015
 *
 */
public abstract class BBTool_ST {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * Must be overridden; the commented body is an example.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		throw new RuntimeException("This method must be overridden.");
		/*
		//Example:
		Timer t=new Timer();
		t.start();
		BBTool_ST bbt=new BBTool_ST(args);
		mb.process(t);
		*/
	}
	
	/**
	 * Constructor.  Parses argument list and sets relevant fields.
	 * Must be called by subclass.
	 * @param args Command line arguments
	 */
	public BBTool_ST(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
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
			}else if(parseArgument(arg, a, b)){
				// do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
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
	
	/**
	 * Must be overridden; this body is just for example.
	 * Parses an argument from the command line.
	 * Assumed to be in "key=value" form, but this is not required.
	 * @param arg The full original argument.
	 * @param a Left hand side, to lower case.
	 * @param b Right hand side, unaltered.
	 * @return true if a matched some keyword.
	 */
	public abstract boolean parseArgument(String arg, String a, String b);

	/**
	 * And example of how to override parseArgument.
	 */
	private boolean parseArgument_EXAMPLE(String arg, String a, String b){
		if(true){throw new RuntimeException("parseArgument() must be overridden.");}
		
		//These are dummy values for demonstration purposes.
		//In real code they should be class fields. 
		int value1;
		boolean value2;
		String value3;
		
		if(a.equals("key1")){
			value1=(int)Tools.parseKMG(b);
			//do anything else necessary here
			return true;
		}else if(a.equals("key2")){
			value2=Tools.parseBoolean(b);
			//do anything else necessary here
			return true;
		}else if(a.equals("key3")){
			value3=b;
			//do anything else necessary here
			return true;
		}
		
		//There was no match to the argument
		return false;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Start the read streams
		startup();
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read streams
		processInner(cris_primary, ros_primary);
		
		//Close the read streams
		shutdown(cris_primary, ros_primary, t);
	}

	/** Create read streams */
	void startup(){
		final ConcurrentReadInputStream cris;
		final ConcurrentReadOutputStream ros;
		
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, false, ffin1, ffin2, qfin1, qfin2);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		final boolean paired=cris.paired();
		outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
		
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
		
		readsProcessed=0;
		basesProcessed=0;
		cris_primary=cris;
		ros_primary=ros;
	}
	
	void shutdown(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final Timer t){
		final boolean paired=cris.paired();
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		errorState|=ReadStats.writeAll(paired);
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads.
	 * This may optionally be overridden. */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					
					boolean keep=processReadPair(r1, r2);
					if(!keep){reads.set(idx, null);}
					
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Process a single read pair.  Must be overidden.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	abstract boolean processReadPair(final Read r1, final Read r2);
	
	/** This is called if the program runs with no parameters.
	 * Should be overidden to print usage information. */
	void printOptions(){
		throw new RuntimeException("printOptions: TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
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
	
	private ConcurrentReadInputStream cris_primary;
	private ConcurrentReadOutputStream ros_primary;
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	protected final FileFormat ffin1;
	protected final FileFormat ffin2;

	protected final FileFormat ffout1;
	protected final FileFormat ffout2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	protected PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	protected boolean overwrite=false;
	protected boolean append=false;
	
}
