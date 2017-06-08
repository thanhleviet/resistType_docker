package jgi;

import java.util.ArrayList;

import align2.ListNum;
import dna.Timer;
import stream.ByteBuilder;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/**
 * Fuses sequences together, with N-padding in between.
 * @author Brian Bushnell
 * @date Jan 20, 2015
 *
 */
public final class FuseSequence extends BBTool_ST {
	
	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		FuseSequence fs=new FuseSequence(args);
		fs.process(t);
	}
	
	public FuseSequence(String[] args){
		super(args);
	}
	
	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#parseArgument(java.lang.String, java.lang.String, java.lang.String)
	 */
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		if(a.equals("pad") || a.equals("npad") || a.equals("ns")){
			npad=Integer.parseInt(b);
			return true;
		}else if(a.equals("q") || a.equals("quality")){
			defaultQuality=Byte.parseByte(b);
			return true;
		}
		return false;
	}
	
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
					
					processReadPair(r1, r2);
					
				}

				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(ros!=null){
			Read r=new Read(bases.toBytes(), quals.toBytes(), 0);
			r.id=(name==null ? "0" : name);
			ArrayList<Read> reads=new ArrayList<Read>(1);
			reads.add(r);
			ros.add(reads, 0);
		}
	}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#processReadPair(stream.Read, stream.Read)
	 */
	@Override
	boolean processReadPair(Read r1, Read r2) {
		if(r1!=null && r1.length()>0){processRead(r1);}
		if(r2!=null && r2.length()>0){processRead(r2);}
		return false;
	}
	
	private void processRead(Read r) {
		if(name==null){name=r.id;}
		if(bases.length>0){
			for(int i=0; i<npad; i++){
				bases.append('N');
				quals.append((byte)0);
			}
		}
		bases.append(r.bases);
		if(r.quality!=null){
			quals.append(r.quality);
		}else{
			for(int i=0, max=r.length(); i<max; i++){
				quals.append(defaultQuality);
			}
		}
	}
	
	int npad=300;
	byte defaultQuality;
	ByteBuilder bases=new ByteBuilder();
	ByteBuilder quals=new ByteBuilder();
	String name;
	
}
