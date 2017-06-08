package kmer;

import java.util.Arrays;

import dna.CoverageArray;


import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;

import align2.IntList;
import align2.LongList;
import align2.Tools;

/**
 * Stores kmers in a long[] and values in an int[][], with a victim cache.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
 */
public abstract class HashArray extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	HashArray(int initialSize, boolean autoResize_, boolean twod){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long)(sizeLimit=(long)(maxLoadFactor*prime));
		array=new long[prime+extra];
		victims=new HashForest(Tools.max(10, initialSize/8), autoResize_, twod);
		Arrays.fill(array, -1);
		autoResize=autoResize_;
		TWOD=twod;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, new int[1]);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Fetched "+Arrays.toString(old));}
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+v+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//	
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			final int[] singleton=new int[1];
//			int[] old=getValues(kmer, singleton);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Before: old="+Arrays.toString(old)+", v="+Arrays.toString(v));}
//			x=set0(kmer, v);
//			if(verbose){System.err.println("After:  old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", get()="+Arrays.toString(getValues(kmer, singleton)));}
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}
//	
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}
	
	@Override
	public final int set(final long kmer, final int[] v){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				if(verbose){System.err.println("A2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				return 0;
			}else if(n==-1){
				if(verbose){System.err.println("B2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
				array[cell]=kmer;
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("B2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				size++;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
		if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
		final int x=victims.set(kmer, v);
		if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
		return x;
	}
	
	@Override
	public final int set(final long kmer, final int v){
		int cell=(int)(kmer%prime);
		
//		assert(TESTMODE);
//		ll.add(kmer);
//		il.add(v);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				if(verbose){System.err.println("A1: Adding "+kmer+", "+v+", "+cell);}
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("A1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				return 0;
			}else if(n==-1){
				if(verbose){System.err.println("B1: Adding "+kmer+", "+v+", "+cell);}
				array[cell]=kmer;
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("B1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				size++;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
		if(verbose){System.err.println("C1: Adding "+kmer+", "+v+", "+cell+
				"; victims.get(kmer)="+Arrays.toString(victims.getValues(kmer, new int[1])));}
		final int x=victims.set(kmer, v);
		if(verbose){System.err.println("C1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1]))+
				"; victims.get(kmer)="+Arrays.toString(victims.getValues(kmer, new int[1])));}
		return x;
	}


//	protected LongList ll=new LongList(); //123
//	protected IntList il=new IntList();
	
	@Override
	public final int setIfNotPresent(long kmer, int value){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				return 0;
			}else if(n==-1){
				array[cell]=kmer;
				insertValue(kmer, value, cell);
				size++;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
//		System.err.println("size="+size+", prime="+prime+", limit="+sizeLimit);
		return victims.setIfNotPresent(kmer, value);
	}
	
	@Override
	public final int getValue(long kmer){
		int cell=findKmer(kmer);
		if(cell==-1){return victims.getValue(kmer);}
		return readCellValue(cell);
	}
	
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		int cell=findKmer(kmer);
		if(cell==-1){return victims.getValues(kmer, singleton);}
		return readCellValues(cell, singleton);
	}
	
	@Override
	public final boolean contains(long kmer){
		int cell=findKmer(kmer);
		if(cell==-1){return victims.contains(kmer);}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	protected abstract void insertValue(final long kmer, final int v, final int cell);
	
	protected abstract void insertValue(final long kmer, final int[] vals, final int cell);

	protected abstract int readCellValue(int cell);
	protected abstract int[] readCellValues(int cell, int[] singleton);
	
	@Override
	final Object get(long kmer){
		throw new RuntimeException("Unimplemented.");
	}
	
	final int findKmer(long kmer){
		int cell=(int)(kmer%prime);
		for(final int max=cell+extra; cell<max; cell++){
			if(array[cell]==kmer){return cell;}
		}
		return -1;
	}
	
	final int findKmerOrEmpty(long kmer){
		int cell=(int)(kmer%prime);
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array[cell];
			if(n==kmer || n==-1){return cell;}
		}
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final boolean canResize() {return true;}
	
	@Override
	final public long size() {return size;}
	
	@Override
	final public int arrayLength() {return array.length;}
	
	@Override
	protected abstract void resize();
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount){
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=-1){
					tsw.print(toText(kmer, readCellValues(i, singleton), k).append('\n'));
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=-1 && readCellValue(i)>=mincount){
					tsw.print(toText(kmer, readCellValue(i), k).append('\n'));
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k, mincount);
		}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount){
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=-1){
					bsw.printlnKmer(kmer, readCellValues(i, singleton), k);
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=-1 && readCellValue(i)>=mincount){
					bsw.printlnKmer(kmer, readCellValue(i), k);
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes(bsw, k, mincount);
		}
		return true;
	}
	
	@Override
	public final void fillHistogram(CoverageArray ca, int max){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=-1){
				int count=Tools.min(readCellValue(i), max);
				ca.increment(count);
			}
		}
		if(victims!=null){
			victims.fillHistogram(ca, max);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long[] array;
	int prime;
	long size=0;
	long sizeLimit;
	final HashForest victims;
	final boolean autoResize;
	public final boolean TWOD;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final static int extra=21;
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE-extra);
	final static float resizeMult=2f; //Resize by a minimum of this much
	final static float minLoadFactor=0.58f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=0.905f; //Reaching this load triggers resizing
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	
	
}
