package stream;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import align2.ListNum;
import align2.Shared;
import align2.Tools;

import fileIO.FileFormat;

public class ConcurrentSolidInputStream extends ConcurrentReadInputStream {
	
	public static void main(String[] args){
		ConcurrentSolidInputStream mates=null;
		if(args.length>2){
			mates=new ConcurrentSolidInputStream(args[2], args[3], 0, null);
		}
		
		ConcurrentSolidInputStream stream=new ConcurrentSolidInputStream(args[0], args[1], 0, mates);
		new Thread(stream).start();
//		ArrayList<Read> list=stream.nextList();
//		System.out.println(list.size());
//		for(int i=0; i<30 && i<list.size(); i++){
//			System.out.println(list.get(i));
//		}
		
		ListNum<Read> ln=stream.nextList();
		long total=0, bad=0;
		while(ln!=null && ln.list!=null && !ln.list.isEmpty()){
			for(int i=0; i<ln.list.size(); i++){
				total++;
//				System.out.println(list.get(i));
				Read r=ln.list.get(i);
				boolean OK=true;
				for(int j=0; j<r.length() && OK; j++){
					if(r.bases[j]=='N'){OK=false;}
				}
				if(!OK){bad++;}
			}
			stream.returnList(ln.id, false);
			ln=stream.nextList();
			if(total>100000){break;}
		}
		
		System.out.println("Total = \t"+total+"\nBad =   \t"+bad+"\t("+String.format("%.3f", bad*100f/total)+"%)");

		stream.shutdown();
		try {
			Thread.sleep(1500);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		}
		for(Thread t : stream.threads){
			System.out.println(t.isAlive());
		}
	}
	
	public ConcurrentSolidInputStream(FileFormat ff1, String qf1, long maxReadsToGenerate, ConcurrentSolidInputStream mateStream_){
		this(ff1.name(), qf1, maxReadsToGenerate, mateStream_);
	}
	
	public ConcurrentSolidInputStream(String fastaName_, String qualName_, long maxReadsToGenerate, 
			ConcurrentSolidInputStream mateStream_){
		fastaName=fastaName_;
		qualName=qualName_;
		producerFasta=new FastaStream(fastaName);
		producerQual=new QualStream(qualName);
		fdepot=new ConcurrentDepot<byte[][]>(BUF_LEN, NUMLISTS_RAW);
		qdepot=new ConcurrentDepot<byte[][]>(BUF_LEN, NUMLISTS_RAW);
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUMLISTS_READ);
		maxReads=maxReadsToGenerate>0 ? maxReadsToGenerate : Long.MAX_VALUE;
		mateStream=mateStream_;
		paired=(mateStream!=null);
	}
	
	public ConcurrentSolidInputStream(FileFormat ff1, String qf1, FileFormat ff2, String qf2, long maxReadsToGenerate){
		this(ff1.name(), qf1, ff2==null ? null : ff2.name(), qf2, maxReadsToGenerate);
	}
	
	public ConcurrentSolidInputStream(String fastaName_, String qualName_, String fastaName2_, String qualName2_, 
			long maxReadsToGenerate){
		fastaName=fastaName_;
		qualName=qualName_;
		producerFasta=new FastaStream(fastaName);
		producerQual=new QualStream(qualName);
		fdepot=new ConcurrentDepot<byte[][]>(BUF_LEN, NUMLISTS_RAW);
		qdepot=new ConcurrentDepot<byte[][]>(BUF_LEN, NUMLISTS_RAW);
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUMLISTS_READ);
		maxReads=maxReadsToGenerate>0 ? maxReadsToGenerate : Long.MAX_VALUE;
		mateStream=(fastaName2_==null ? null : new ConcurrentSolidInputStream(fastaName2_, qualName2_, maxReadsToGenerate, null));
		paired=(mateStream!=null);
	}
	
	/** If running in paired-end mode, attaches mated reads to each other. */
	private ArrayList<Read> makeReadList() {
		ArrayList<Read> rlist=makeReadList2();
		if(mateStream!=null){
			ListNum<Read> matesln=mateStream.nextList();
			ArrayList<Read> mates=matesln.list;
			if(rlist!=null && mates!=null){
				int max=Tools.min(rlist.size(), mates.size());
				for(int i=0; i<max; i++){
					Read a=rlist.get(i);
					Read b=mates.get(i);
					a.mate=b;
					b.mate=a;
					b.setPairnum(1);
				}
				mates.clear();
				mateStream.returnList(matesln.id, false);
			}
		}
		return rlist;
	}
	
	private ArrayList<Read> makeReadList2() {
		
		ArrayList<byte[][]> flist=null;
		ArrayList<byte[][]> qlist=null;
		ArrayList<Read> rlist=null;
//		System.out.println("Making list");
		while(flist==null){
			try {
				flist=fdepot.full.take();
//				System.out.println("Got flist of size "+flist.size());
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				if(shutdown){return null;}
				//e.printStackTrace();
			}
		}
		while(qlist==null){
			try {
				qlist=qdepot.full.take();
//				System.out.println("Got qlist of size "+flist.size());
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				if(shutdown){return null;}
				//e.printStackTrace();
			}
		}
		while(rlist==null){
			try {
				rlist=depot.empty.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				if(shutdown){return null;}
				//e.printStackTrace();
			}
		}
		
		assert(flist.size()==qlist.size() || 
				(maxReads<Long.MAX_VALUE && (flist.isEmpty() || qlist.isEmpty()))) : flist.size()+", "+qlist.size()+", "+shutdown;
		if(flist.isEmpty() || qlist.isEmpty()){
//			System.out.println("Returning early... ");
			fdepot.full.add(flist);
			qdepot.full.add(qlist);
			return rlist;
		}
		
		for(int i=0; i<flist.size() && generated<maxReads; i++){
			byte[][] fasta=flist.get(i);
			byte[][] quals=qlist.get(i);
			Read r=new Read(fasta, quals, true, generated);
			readsIn++;
			basesIn+=r.length();
			if(randy==null || randy.nextFloat()<samplerate){
				rlist.add(r);
			}
			generated++;
			if(generated>1 && ((generated%1000000)==0) && mateStream==null){System.err.println("Generated read #"+generated);}
		}
		flist.clear();
		qlist.clear();
		fdepot.empty.add(flist);
		qdepot.empty.add(qlist);
		
		return rlist;
	}
	
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> rlist=null;
		
		while(rlist==null){
			try {
				//System.err.println((mateStream==null ? 2 : 1)+" Attempting take; depot.size = "+rdepot.full.size()+", "+rdepot.empty.size());
				rlist=depot.full.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				//e.printStackTrace();
			}
		}
		//System.err.println((mateStream==null ? 2 : 1)+" Took "+rlist.size());
		ListNum<Read> ln=new ListNum<Read>(rlist, listnum);
		listnum++;
		return ln;
	}
	
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("crisS:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("crisS:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));
		}
	}
	
	@Override
	public void run() {
		FThread fthread=new FThread();
		QThread qthread=new QThread();
		threads=new Thread[] {new Thread(fthread), new Thread(qthread), Thread.currentThread()};
		threads[0].start();
		threads[1].start();

		if(mateStream!=null){new Thread(mateStream).start();}
		
		ArrayList<Read> list=makeReadList();
		while(list!=null && !list.isEmpty() && !shutdown){
			//System.err.println((mateStream==null ? 2 : 1)+" Adding list to rdepot size "+rdepot.full.size()+"/"+rdepot.full.remainingCapacity());
			depot.full.add(list);
			list=makeReadList();
		}
		
		//System.err.println((mateStream==null ? 2 : 1)+" Exiting main loop.");
		if(generated>=maxReads){shutdown();}
		//System.err.println((mateStream==null ? 2 : 1)+" Shutdown complete.");
		
		if(list!=null){
			//System.err.println((mateStream==null ? 2 : 1)+" Attempting to add current list to rdepot.full: "+rdepot.full.size());
			list.clear();
			depot.full.add(list);
		}else{
			assert(shutdown) : "Null read list encountered for unknown reason.";
//			System.err.println("Null read list encountered.");
//			shutdown();
		}
		
		//Add poison pills
		//System.err.println((mateStream==null ? 2 : 1)+" Attempting to add poison list to rdepot.full: "+rdepot.full.size());
		depot.full.add(new ArrayList<Read>());
		
		for(int i=1; i<depot.bufferCount; i++){
			list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					//e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){depot.full.add(list);}else{break;}
		}
		//End thread
		
		while(!depot.empty.isEmpty()){
			depot.full.add(depot.empty.poll());
		}
//		System.err.println(rdepot.full.size()+", "+rdepot.empty.size());

		
		if(mateStream!=null){mateStream.shutdown();}
	}
	
	private class FThread implements Runnable{

		@Override
		public void run() {
			while(!shutdown && producerFasta.hasMore()){
				ArrayList<byte[][]> list=null;
				while(list==null){
					try {
						list=fdepot.empty.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						//e.printStackTrace();
						if(shutdown){return;}
					}
				}
				for(int i=0; i<fdepot.bufferSize; i++){
					byte[][] r=producerFasta.next();
					if(r==null){break;}
					list.add(r);
				}
				fdepot.full.add(list);
				fgenerated+=list.size();
				if(fgenerated>=maxReads){break;}
			}
//			System.err.println("Shutting down FThread.");
			
			//Add poison pills
			for(int i=1; i<fdepot.bufferCount; i++){
				ArrayList<byte[][]> list=null;
				while(list==null){
					try {
						list=fdepot.empty.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						//e.printStackTrace();
						if(shutdown){return;}
					}
				}
				fdepot.full.add(list);
			}
//			System.err.println("Done shutting down FThread.");
			
			//End thread
			producerFasta.close();
		}

		private long fgenerated=0;
	}
	
	private class QThread implements Runnable{

		@Override
		public void run() {
			while(!shutdown && producerQual.hasMore()){
				ArrayList<byte[][]> list=null;
				while(list==null){
					try {
						list=qdepot.empty.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						//e.printStackTrace();
						if(shutdown){return;}
					}
				}
				for(int i=0; i<qdepot.bufferSize; i++){
					byte[][] r=producerQual.next();
					if(r==null){break;}
					list.add(r);
				}
				qdepot.full.add(list);
				qgenerated+=list.size();
				if(qgenerated>=maxReads){break;}
			}
//			System.err.println("Shutting down QThread.");
			
			//Add poison pills
			for(int i=1; i<qdepot.bufferCount; i++){
				ArrayList<byte[][]> list=null;
				while(list==null){
					try {
						list=qdepot.empty.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						//e.printStackTrace();
						if(shutdown){return;}
					}
				}
				qdepot.full.add(list);
			}
			producerQual.close();
			//End thread
//			System.err.println("Done shutting down QThread.");
		}
		
		private long qgenerated=0;
	}

	public void shutdown(){
		synchronized(ShutdownKey){
			if(shutdown){return;}
			//System.err.println("Shutting down SCRIS.");
			shutdown=true;
			ShutdownKey[0]=true;
			//System.err.println("A");
			threads[0].interrupt();
			//System.err.println("B");
			threads[1].interrupt();
			//System.err.println("C");
			threads[2].interrupt();
		}
		if(mateStream!=null){
			mateStream.shutdown();
		}
	}

	@Override
	public synchronized void restart() {
		shutdown=false;
		generated=0;	
		basesIn=0;
		readsIn=0;
		producerFasta=new FastaStream(fastaName);
		producerQual=new QualStream(qualName);
		fdepot=new ConcurrentDepot<byte[][]>(BUF_LEN, NUMLISTS_RAW);
		qdepot=new ConcurrentDepot<byte[][]>(BUF_LEN, NUMLISTS_RAW);
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUMLISTS_READ);
	}

	@Override
	public synchronized void close() {
		producerFasta.close();
		producerQual.close();
	}
	
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		if(rate>=1f){
			randy=null;
		}else if(seed>-1){
			randy=new java.util.Random(seed);
		}else{
			randy=new java.util.Random();
		}
	}
	
	public long basesIn(){return basesIn;}
	public long readsIn(){return readsIn;}
	
	@Override
	public boolean paired() {
		assert(paired==(mateStream!=null));
		return paired;
	}
	
	@Override
	public boolean verbose(){return verbose;}
	
	private float samplerate=1f;
	private java.util.Random randy=null;
	
	public Object[] producers(){return new Object[] {producerFasta, producerQual};}
	
	@Override
	public boolean errorState(){return errorState;}
	/** TODO */
	private boolean errorState=false;
	
	private boolean shutdown=false;

	private Thread[] threads;
	
	private long basesIn=0;
	private long readsIn=0;
	
	public final long maxReads;
	private long generated=0;
	private long listnum=0;

	public final String fastaName;
	public final String qualName;
	
	public FastaStream producerFasta;
	public QualStream producerQual;
	private ConcurrentDepot<byte[][]> fdepot;
	private ConcurrentDepot<byte[][]> qdepot;
	private ConcurrentDepot<Read> depot;

	public static int NUMLISTS_RAW=7;
	public static int NUMLISTS_READ=25;
	public static boolean verbose=false;
	
	private final ConcurrentSolidInputStream mateStream;
	private final boolean paired;

//	private final Object ShutdownKey="ShutdownKey"+hashCode();
	private final boolean[] ShutdownKey=new boolean[] {false};
	
}
