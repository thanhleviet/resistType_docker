package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ByteBuilder;
import align2.LongList;
import align2.ReadStats;
import align2.Tools;

import dna.Parser;
import dna.Timer;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date Oct 28, 2014
 *
 */
public class CallPeaks {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		CallPeaks cp=new CallPeaks(args);
		cp.process(t);
	}
	
	public CallPeaks(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		if(printClass){outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");}
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(Parser.isJavaFlag(arg)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				in=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("minheight") || a.equals("h")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume") || a.equals("v")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth") || a.equals("w")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak") || a.equals("minp")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak") || a.equals("maxp")){
				maxPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(in==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		
		if(out==null){out="stdout.txt";}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TEXT, null, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		LongList hist=loadHistogram(ffin);
		ArrayList<Peak> peaks=callPeaks(hist);
		hist=null;
		printPeaks(peaks);
		t.stop();
		System.err.println("\nFound "+peaks.size()+" peaks in "+t);
		
		peaks=null;
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	public static boolean printPeaks(long[] array, String fname, boolean ow, long minHeight, long minVolume, int minWidth, int minPeak, int maxPeak, int maxPeakCount){
		CallPeaks cp=new CallPeaks(new String[] {"out="+fname, "ow="+ow,
				"minheight="+minHeight, "minvolume="+minVolume, "minwidth="+minWidth, "minpeak="+minPeak, "maxpeak="+maxPeak, "maxpeaks="+maxPeakCount});
		ArrayList<Peak> peaks=cp.callPeaks(array, array.length);
		cp.printPeaks(peaks);
		return cp.errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static LongList loadHistogram(FileFormat ff){
		LongList list=new LongList(8000);
		TextFile tf=new TextFile(ff);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("#")){
				//ignore
			}else{
				String[] split=line.split("\\s+");
				if(split.length==1){
					list.add(Long.parseLong(split[0]));
				}else{
					list.set(Integer.parseInt(split[0]), Long.parseLong(split[1]));
				}
			}
		}
		boolean errorState_=tf.close();
		assert(!errorState_) : "Encountered an error when reading "+ff.name()+".\n" +
				"To skip this error message, run with the '-da' flag.";
		
		return list;
	}
	
	private static ArrayList<Peak> condense(ArrayList<Peak> in, int maxCount){
		if(in==null || in.isEmpty()){return in;}
		maxCount=Tools.max(Tools.min(in.size(), maxCount), 1);
		ArrayList<Peak> out=new ArrayList<Peak>(Tools.min(maxCount, in.size()));
		long[] counts=new long[in.size()];
		for(int i=0; i<in.size(); i++){
			Peak p=in.get(i);
			counts[i]=(callByRawCount ? p.centerHeight2() : p.centerHeight);
		}
		Arrays.sort(counts);
		long limit=counts[counts.length-maxCount];
		for(Peak p : in){
			if((callByRawCount ? p.centerHeight2() : p.centerHeight)>=limit){out.add(p);}
		}
		for(Peak p : in){
			if((callByRawCount ? p.centerHeight2() : p.centerHeight)<limit){
				Peak p2=out.get(0);
				for(Peak temp : out){
					if(Tools.absdif(p.center, temp.center)<Tools.absdif(p.center, p2.center)){
						p2=temp;
					}
				}
				p2.absorb(p);
			}
		}
		return out;
	}
	
	private void printPeaks(ArrayList<Peak> peaks){
		if(ffout==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		bsw.println("#start\tcenter\tstop\tmax\tvolume");
		ByteBuilder bb=new ByteBuilder(200);
		for(Peak p : peaks){
			p.toBytes(bb);
			bsw.println(bb);
			bb.setLength(0);
		}
		errorState|=bsw.poisonAndWait();
	}
	
	public ArrayList<Peak> callPeaks(LongList list){
		return callPeaks(list.array, list.size);
	}
	
	public ArrayList<Peak> callPeaks(long[] array, int length){
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		
		int dip0=-1;
		for(int i=1; i<length; i++){
			if(array[i-1]<array[i]){
				dip0=i-1;
				break;
			}
		}
		if(dip0<0){return peaks;}
//		assert(false) : dip0+", "+array[dip0);
		
		final int UP=0, DOWN=1;
		int mode=UP;
		int start=dip0, center=-1;
		long prev=array[dip0];
		long sum=prev;
		long sum2=prev*dip0;
		for(int i=dip0+1; i<length; i++){
			final long x=array[i];
			
//			if(i<16){System.err.println("i="+i+", x="+x+", mode="+mode+", center="+center+", start="+start+", dip0="+dip0);}
			
			if(mode==UP){
				if(x<prev){
					mode=DOWN;
					center=i-1;
				}
			}else{
				if(x>prev){
					mode=UP;
					int stop=i-1;
					long max=array[center];
					if(center>=minPeak && center<=maxPeak && max>=minHeight && (stop-start)>=minWidth && sum>=minVolume){
						for(int j=center-1; j>=0; j--){//find middle of mesas
							if(array[j]!=max){
								center=(center+j+2)/2;
								break;
							}
						}
						{
							long valley=array[stop];
							for(int j=stop; j>=0; j--){//find middle of valleys
								if(array[j]!=valley){
									if(valley==0){stop=j+1;}
									else{stop=(stop+j+2)/2;}
									break;
								}
							}
						}
						
						Peak p=new Peak(center, start, stop, max, array[start], array[stop], sum, sum2);
						peaks.add(p);
					}else{
//						Peak p=new Peak(center, start, stop, max, sum);
//						System.err.println("*"+p);
					}
					start=stop;
					stop=-1;
					sum=sum2=0;
					center=-1;
					if(i>maxPeak){break;}
					while(i<array.length && array[i]==0){i++;}//Skip zero regions
				}
			}
			
			sum+=x;
			sum2+=(x*i);
			prev=x;
		}
		
		if(mode==DOWN){
			int stop=length;
			long max=array[center];
			for(int j=center-1; j>=0; j--){//find middle of mesas
				if(array[j]!=max){
					center=(center+j+2)/2;
					break;
				}
			}
			{
				long valley=array[stop-1];
				for(int j=stop-1; j>=0; j--){//find middle of valleys
					if(array[j]!=valley){
						if(valley==0){stop=j+1;}
						else{stop=(stop+j+2)/2;}
						break;
					}
				}
			}
			if(center>=minPeak && center<=maxPeak && max>=minHeight && (stop-start)>=minWidth && sum>=minVolume){
				Peak p=new Peak(center, start, stop, max, array[start], array[Tools.min(stop, length-1)], sum, sum2);
				peaks.add(p);
			}else{
//				Peak p=new Peak(center, start, stop, max, sum);
//				System.err.println("*"+p);
			}
		}
		
		if(maxPeakCount<peaks.size()){
			peaks=condense(peaks, maxPeakCount);
		}
		return peaks;
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		throw new RuntimeException("printOptions: TODO");
//		outstream.println("Syntax:\n");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class Peak{
		
		Peak(int center_, int start_, int stop_, long centerHeight_, long startHeight_, long stopHeight_, long volume_, long volume2_){
			
			center=center_;
			start=start_;
			stop=stop_;
			
			centerHeight=centerHeight_;
			startHeight=startHeight_;
			stopHeight=stopHeight_;
			volume=volume_;
			volume2=volume2_;
			
			assert(center>=0) : this;
			assert(start<center) : this;
			assert(stop>center) : this;
		}
		
		/**
		 * @param p
		 */
		public void absorb(Peak p) {
			
			if(center>p.center){
				assert(p.stop<stop) : "\n"+this+"\n"+p+"\n";
				if(start>p.start){
					start=p.start;
					startHeight=p.startHeight;
				}
			}else{
				assert(p.start>start) : "\n"+this+"\n"+p+"\n";
				if(stop<p.stop){
					stop=p.stop;
					stopHeight=p.stopHeight;
				}
			}
			
			long c1=callByRawCount ? centerHeight2() : centerHeight;
			long c2=callByRawCount ? p.centerHeight2() : p.centerHeight;
			
			if(c1<c2){
				center=p.center;
				centerHeight=p.centerHeight;
			}
			
			volume+=p.volume;
			volume2+=p.volume2;
			
		}

		int width(){return stop-start;}
		
		@Override
		public String toString(){
			return start+"\t"+center+"\t"+stop+"\t"+centerHeight+"\t"+volume;
		}
		
		public ByteBuilder toBytes(ByteBuilder bb){
			if(bb==null){bb=new ByteBuilder();}
			bb.append(start);
			bb.append('\t');
			bb.append(center);
			bb.append('\t');
			bb.append(stop);
			bb.append('\t');
			bb.append(centerHeight);
			bb.append('\t');
			bb.append(volume);
			bb.append('\t');
			return bb;
		}

		/** Inclusive */
		public int start;
		public int center;
		/** Exclusive */
		public int stop;

		//Unique counts
		public long startHeight;
		public long centerHeight;
		public long stopHeight;
		public long volume;

		public long volume2;

		//Raw counts
		public long startHeight2(){return startHeight*start;}
		public long centerHeight2(){return centerHeight*center;}
		public long stopHeight2(){return stopHeight*stop;}
		
		
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long minHeight=2;
	private long minVolume=2;
	private int minWidth=2;
	private int minPeak=2;
	private int maxPeak=Integer.MAX_VALUE;
	private int maxPeakCount=8;
	
	private String in;
	private String out;
	
	private final FileFormat ffin;
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	public static boolean printClass=true;
	
	public static boolean callByRawCount=true;
	
}
