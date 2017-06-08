package kmer;

import dna.CoverageArray;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class HashBuffer extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashBuffer(AbstractKmerTable[] tables_, int buflen_){
		tables=tables_;
		buflen=buflen_;
		assert(buflen<Short.MAX_VALUE);
		ways=tables.length;
		sizes=new short[ways];
		buffers=new long[ways][buflen];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int incrementAndReturnNumCreated(long kmer) {
		final int way=(int)(kmer%ways);
		long[] buffer=buffers[way];
		buffer[sizes[way]]=kmer;
		sizes[way]++;
		if(sizes[way]>=buflen){
			return incrementBuffer(way);
		}
		return 0;
	}
	
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=incrementBuffer(i);}
		return added;
	}
	
	@Override
	public int set(long kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int set(long kmer, int[] vals) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int setIfNotPresent(long kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int getValue(long kmer) {
		final int way=(int)(kmer%ways);
		return tables[way].getValue(kmer);
	}
	
	@Override
	public int[] getValues(long kmer, int[] singleton){
		final int way=(int)(kmer%ways);
		return tables[way].getValues(kmer, singleton);
	}
	
	@Override
	public boolean contains(long kmer) {
		final int way=(int)(kmer%ways);
		return tables[way].contains(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	Object get(long kmer) {
		final int way=(int)(kmer%ways);
		return tables[way].get(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private int incrementBuffer(final int way){
		final int size=sizes[way];
		if(size<1){return 0;}
		sizes[way]=0;
		final long[] buffer=buffers[way];
		int added=0;
		final AbstractKmerTable table=tables[way];
		synchronized(table){
			for(int i=0; i<size; i++){
				final long kmer=buffer[i];
				added+=table.incrementAndReturnNumCreated(kmer);
			}
		}
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final boolean canResize() {return false;}
	
	@Override
	public final boolean canRebalance() {return false;}
	
	@Deprecated
	@Override
	public long size() {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Unimplemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount){
		for(AbstractKmerTable table : tables){
			dumpKmersAsText(tsw, k, mincount);
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount){
		for(AbstractKmerTable table : tables){
			dumpKmersAsBytes(bsw, k, mincount);
		}
		return true;
	}
	
	@Override
	public void fillHistogram(CoverageArray ca, int max){
		for(AbstractKmerTable table : tables){
			fillHistogram(ca, max);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int increment(long kmer) {
		throw new RuntimeException("Unsupported");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final AbstractKmerTable[] tables;
	private final int buflen;
	private final int ways;	
	private final short[] sizes;
	private final long[][] buffers;

}
