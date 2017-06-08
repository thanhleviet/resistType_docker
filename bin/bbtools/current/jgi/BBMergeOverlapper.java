package jgi;

import java.util.Arrays;
import java.io.File;

import kmer.KCountArray;
import stream.Read;
import align2.QualityTools;
import align2.Tools;
import dna.AminoAcid;
import align2.Shared;

/**
 * @author Brian Bushnell
 * @date Apr 15, 2014
 *
 */
public final class BBMergeOverlapper {

	static {
		if(Shared.USE_JNI){
			String name = "bbtoolsjni";
			try {
				System.loadLibrary(name);
			} catch (UnsatisfiedLinkError e1) {
				// System.loadLibrary() does not work with MPI.
				// Need to use System.load() with an explicit full
				// path to the native library file for the MPI case.
				boolean success = false;
				String libpath=System.getProperty("java.library.path");
				libpath = libpath.replace("-Djava.library.path=","");
				String[] libpathEntries = libpath.split(File.pathSeparator);
				for(int i = 0; i < libpathEntries.length; i++) {
					if(success) break;
					String lib = libpathEntries[i]+"/"+System.mapLibraryName(name);
					try {
						System.load(lib);
						success = true;
					} catch (UnsatisfiedLinkError e2) {
						success = false;
						if((i+1) >= libpathEntries.length) {
							throw new RuntimeException("\n\n*****   Native library can not be found in java.library.path.   *****\n");
						}
					}
				}
			}
		}
	}
	
	private static native final int mateByOverlapJNI(byte[] a_bases, byte[] b_bases, byte[] a_quality, byte[] b_quality,
			int[] rvector, int minOverlap0, int minOverlap, int margin,
			int maxMismatches0, int maxMismatches, int minq, int a_insertSizeMapped);
	
	protected static final int mateByOverlap(Read a, Read b, int[] rvector, int minOverlap0, final int minOverlap, int margin, final int maxMismatches0, final int maxMismatches, final int minq) {
		if(Shared.USE_JNI){
			return mateByOverlapJNI(a.bases,b.bases,a.quality,b.quality,rvector,minOverlap0,minOverlap,margin,maxMismatches0,maxMismatches,minq,a.insertSizeMapped(IGNORE_MAPPING_STRAND));
		}else{
//			return mateByOverlapJava(a, b, rvector, minOverlap0, minOverlap, margin, maxMismatches0, maxMismatches, minq);
			return mateByOverlapJava_unrolled(a, b, rvector, minOverlap0, minOverlap, margin, maxMismatches0, maxMismatches, minq);
		}
	}
	
	protected static final int mateByOverlapJava(Read a, Read b, int[] rvector, int minOverlap0, final int minOverlap, int margin,
			final int maxMismatches0, final int maxMismatches, final int minq) {
		minOverlap0=Tools.min(Tools.max(1, minOverlap0), minOverlap);
		assert(maxMismatches<=maxMismatches0);
		margin=Tools.max(margin, 0);
		assert(maxMismatches>=margin);
		if(rvector==null){rvector=new int[5];}
		if(USE_MAPPING){
			rvector[0]=100;
			rvector[1]=20;
			rvector[2]=0;
			//rvector[3]=20; Unused
			rvector[4]=0;
			return a.insertSizeMapped(IGNORE_MAPPING_STRAND);
		}
		
//		assert(false) : minOverlap0+", "+minOverlap+", "+margin+", "+maxMismatches0+", "+maxMismatches+", "+minq+", "+MIN_OVERLAP_INSERT;
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		
		int bestOverlap=-1;
//		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=maxMismatches0;
//		float bestExpected=0;
		
		boolean ambig=false;
		final int maxOverlap=abases.length+bbases.length-Tools.max(minOverlap, MIN_OVERLAP_INSERT);
//		assert(false) : minOverlap+", "+maxOverlap;
//		System.err.print("\nm");
		
		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
//			System.err.print("\nn");
//			verbose=(insert==174);
			if(verbose){System.err.println("\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(abases.length+bbases.length-overlap));}
			
			
			int good=0, bad=0;
//			float expected=0;
			
			int istart=(overlap<=abases.length ? 0 : overlap-abases.length);
			int jstart=(overlap<=abases.length ? abases.length-overlap : 0);
//			System.err.print("o");
			
			for(int i=istart, j=jstart, badlim=bestBad+margin; i<overlap && i<bbases.length && j<abases.length && bad<=badlim; i++, j++){
				assert(j>=0 && j<=abases.length && i>=0 && i<=bbases.length) : "\njstart="+jstart+", j="+j+
					", istart="+istart+", i="+i+" \n"+"overlap="+overlap+", a.length="+a.length()+
					", b.length="+b.length()+", bad="+bad+", badlim="+badlim+", good="+good;
				final byte ca=abases[j], cb=bbases[i];
				final byte qa=(aqual==null ? 20 : aqual[j]), qb=(bqual==null ? 20 : bqual[i]);
				byte q=QualityTools.qualsToPhred(qa, qb);
//				if(ca=='N' || cb=='N' || (aqual!=null && aqual[j]<minq) || (bqual!=null && bqual[i]<minq)){
				if(ca=='N' || cb=='N' || q<minq){
					//do nothing
				}else{
					assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
					if(ca==cb){good++;}
					else{bad++;}
//					expected+=QualityTools.qualsToProbError(qa, qb);
				}
			}
//			System.err.print("p");
			
			if(verbose){System.err.println("overlap="+overlap+", bad="+bad+", good="+good+", bestOverlap="+
					bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig);}

//			System.err.print("a");
			if(bad*2<good){
				if(good>minOverlap){//Candidate
					if(bad<=bestBad){

						//					System.err.print("b");
						if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
//							ambig=(bestBad-bad<margin);
							if(bestBad-bad<margin){ambig=true;}
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
//							bestExpected=expected;
						}else if(bad==bestBad){
							ambig=true;
						}

						//					System.err.print("c");
						if(ambig && bestBad<margin){
							//						System.err.print("d");
							rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
							rvector[1]=bestGood;
							rvector[2]=bestBad;
							//						rvector[3]=bestSum;
							rvector[4]=(ambig ? 1 : 0);
							//						System.err.print("e");
							return -1;
						}

						//					System.err.print("f");
					}
				}else if(bad<margin){
					//				System.err.print("g");
					ambig=true;
					rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
					rvector[1]=bestGood;
					rvector[2]=bestBad;
					//				rvector[3]=bestSum;
					rvector[4]=(ambig ? 1 : 0);
					return -1;
				}
			}
//			System.err.print("h");
			
		}
//		System.err.println("i");
		
//		if(bestExpected<bestBad+.5f){ambig=true;}
//		if((bestExpected+.04f)*6f<bestBad){ambig=true;}
		
		if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}
		
		rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
		rvector[1]=bestGood;
		rvector[2]=bestBad;
//		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
//		assert(false) : Arrays.toString(rvector)+"\n"+new String(a.bases)+"\n"+new String(b.bases);
		
		return (bestOverlap<0 ? -1 : abases.length+bbases.length-bestOverlap);
	}
	
	protected static final int mateByOverlapJava_unrolled(Read a, Read b, int[] rvector, int minOverlap0, final int minOverlap, int margin,
			final int maxMismatches0, final int maxMismatches, final int minq) {
		minOverlap0=Tools.min(Tools.max(1, minOverlap0), minOverlap);
		assert(maxMismatches<=maxMismatches0);
		margin=Tools.max(margin, 0);
		assert(maxMismatches>=margin);
		if(rvector==null){rvector=new int[5];}
		if(USE_MAPPING){
			rvector[0]=100;
			rvector[1]=20;
			rvector[2]=0;
			//rvector[3]=20; Unused
			rvector[4]=0;
			return a.insertSizeMapped(IGNORE_MAPPING_STRAND);
		}
		
//		assert(false) : minOverlap0+", "+minOverlap+", "+margin+", "+maxMismatches0+", "+maxMismatches+", "+minq+", "+MIN_OVERLAP_INSERT;
		
		final byte[] abases=a.bases, bbases=b.bases;
		final int alen=abases.length, blen=bbases.length;
		final byte[] aqual=(a.quality==null ? alen<FAKE_LEN ? fakeQual : QualityTools.fakeQuality(FAKE_Q, alen) : a.quality);
		final byte[] bqual=(b.quality==null ? blen<FAKE_LEN ? fakeQual : QualityTools.fakeQuality(FAKE_Q, blen) : b.quality);
		
		int bestOverlap=-1;
//		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=maxMismatches0;
//		float bestExpected=0;
		
		boolean ambig=false;
		final int maxOverlap=alen+blen-Tools.max(minOverlap, MIN_OVERLAP_INSERT);
//		assert(false) : minOverlap+", "+maxOverlap;
//		System.err.print("\nm");
		
		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
//			System.err.print("\nn");
//			verbose=(insert==174);
			if(verbose){System.err.println("\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(alen+blen-overlap));}
			
			
			int good=0, bad=0;
//			float expected=0;
			
			int istart=(overlap<=alen ? 0 : overlap-alen);
			int jstart=(overlap<=alen ? alen-overlap : 0);
//			System.err.print("o");
			
			{
				final int iters=Tools.min(overlap-istart, blen-istart, alen-jstart);
				final int maxi=istart+iters;
				final int badlim=bestBad+margin;
				final int iter=0; //For unrolling
				
//				for(int preloops=iters&1; iter<preloops && bad<=badlim; iter++){
//					final int i=istart+iter;
//					final int j=jstart+iter;
//					assert(j>=0 && j<=alen && i>=0 && i<=blen) : "\njstart="+jstart+", j="+j+
//						", istart="+istart+", i="+i+" \n"+"overlap="+overlap+", a.length="+alen+
//						", b.length="+blen+", bad="+bad+", badlim="+badlim+", good="+good;
//					final byte ca=abases[j], cb=bbases[i];
//					final byte q=QualityTools.qualsToPhred_Capped(aqual[j], bqual[i]);
//
//					if(ca=='N' || cb=='N' || q<minq){
//						//do nothing
//					}else{
//						assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
//						if(ca==cb){good++;}
//						else{bad++;}
//					}
//				}
//				
//				for(int i1=istart+iter, j1=jstart+iter, i2=istart+iter+1, j2=jstart+iter+1; i1<maxi && bad<=badlim; i1+=2, j1+=2, i2+=2, j2+=2){
//					assert(j1>=0 && j1<=alen && i1>=0 && i1<=blen) : "\njstart="+jstart+", j="+j1+
//						", istart="+istart+", i="+i1+" \n"+"overlap="+overlap+", a.length="+alen+
//						", b.length="+blen+", bad="+bad+", badlim="+badlim+", good="+good;
//					assert(j1<alen) : j1+", "+alen+", "+i1+", "+blen+", "+iter+", "+iters+", "+istart+", "+jstart;
//					assert(i1<blen) : j1+", "+alen+", "+i1+", "+blen+", "+iter+", "+iters+", "+istart+", "+jstart;
//					assert(j2<alen) : j2+", "+alen+", "+i2+", "+blen+", "+iter+", "+iters+", "+istart+", "+jstart;
//					assert(i2<blen) : j2+", "+alen+", "+i2+", "+blen+", "+iter+", "+iters+", "+istart+", "+jstart;
//					final byte ca1=abases[j1], cb1=bbases[i1], ca2=abases[j2], cb2=bbases[i2];
//					final byte q1=QualityTools.qualsToPhred_Capped(aqual[j1], bqual[i1]), q2=QualityTools.qualsToPhred_Capped(aqual[j2], bqual[i2]);
//					
//					if(ca1=='N' || cb1=='N' || q1<minq){//do nothing
//					}else if(ca1==cb1){good++;}
//					else{bad++;}
//					
//					if(ca2=='N' || cb2=='N' || q2<minq/* || bad>=badlim*/){//do nothing
//					}else if(ca2==cb2){good++;}
//					else{bad++;}
//				}
				
				for(int i1=istart+iter, j1=jstart+iter; i1<maxi && bad<=badlim; i1++, j1++){
					assert(j1>=0 && j1<=alen && i1>=0 && i1<=blen) : "\njstart="+jstart+", j="+j1+
						", istart="+istart+", i="+i1+" \n"+"overlap="+overlap+", a.length="+alen+
						", b.length="+blen+", bad="+bad+", badlim="+badlim+", good="+good;
					final byte ca1=abases[j1], cb1=bbases[i1];
					final byte q1=QualityTools.qualsToPhred(aqual[j1], bqual[i1]);
					
//					System.err.println(i1+","+j1);
					
					if(ca1=='N' || cb1=='N' || q1<minq){//do nothing
					}else if(ca1==cb1){good++;}
					else{bad++;}
				}
				
				if(verbose){
					final int overlapLen=(maxi-istart);
//					System.err.println(istart+", "+maxi+", "+overlapLen);
//					System.err.println(jstart+", "+(jstart+overlapLen)+", "+overlapLen);
					System.err.println("overlapLen="+overlapLen+"; coordinates ("+jstart+"-"+(jstart+overlapLen)+"), ("+istart+"-"+maxi+")");
					System.err.println(new String(abases, jstart, overlapLen));
					System.err.println(new String(bbases, istart, overlapLen));
				}
				
				if(verbose){System.err.println("overlap="+overlap+", bad="+bad+", good="+good+", badlim="+badlim+", bestOverlap="+
						bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig+", mino="+minOverlap+", mino0="+minOverlap0+
						", margin="+margin+", maxMismatches="+maxMismatches);}	
			}
			
//			System.err.print("p");
			
//			System.err.println(overlap+", "+bestOverlap+", "+bestGood+", "+bestBad+", "+ambig);
			
//			System.err.print("a");
			if(bad*2<good){
				if(verbose){System.err.print("a");}
				if(good>minOverlap){//Candidate
					if(verbose){System.err.print("b");}
					if(bad<=bestBad){

						if(verbose){System.err.print("c");}
						if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
							if(verbose){System.err.print("d");}
//							ambig=(bestBad-bad<margin);
							if(bestBad-bad<margin){ambig=true;}
//							else if(bestBad-bad>margin && good>bestGood){ambig=false;} //Increases FP rate with no real advantage
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
//							bestExpected=expected;
						}else if(bad==bestBad){
							if(verbose){System.err.print("e");}
							ambig=true;
						}
						
						//					System.err.print("c");
						if(ambig && bestBad<margin){
							if(verbose){System.err.print("f");}
							rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
							rvector[1]=bestGood;
							rvector[2]=bestBad;
							//						rvector[3]=bestSum;
							rvector[4]=(ambig ? 1 : 0);
							//						System.err.print("e");
							return -1;
						}

						//					System.err.print("f");
					}
				}else if(bad<margin){
					if(verbose){System.err.print("g");}
					ambig=true;
					rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
					rvector[1]=bestGood;
					rvector[2]=bestBad;
					//				rvector[3]=bestSum;
					rvector[4]=(ambig ? 1 : 0);
					return -1;
				}else{
					if(verbose){System.err.print("h");}
				}
			}
//			System.err.print("h");
			
		}
//		System.err.println("i");
		
//		if(bestExpected<bestBad+.5f){ambig=true;}
//		if((bestExpected+.04f)*6f<bestBad){ambig=true;}
		
		if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}
		
		rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
		rvector[1]=bestGood;
		rvector[2]=bestBad;
//		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
//		assert(false) : Arrays.toString(rvector)+"\n"+new String(a.bases)+"\n"+new String(b.bases);
		
		if(verbose){System.err.println("bestOverlap="+
				bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig+", mino="+minOverlap+", mino0="+minOverlap0+
				", margin="+margin+", maxMismatches="+maxMismatches);}	
		
		return (bestOverlap<0 ? -1 : alen+blen-bestOverlap);
	}
	
	
	/**
	 * TODO Use this
	 * @param a
	 * @param b
	 * @param overlap
	 * @return
	 */
	protected static final float expectedMismatches(Read a, Read b, int overlap) {
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		int good=0, bad=0;
		float expected=0;
		float goodWeighted=0;
		float badWeighted=0;
		
		if(aqual==null || bqual==null){return (overlap+0)/16;}

		int istart=(overlap<=abases.length ? 0 : overlap-abases.length);
		int jstart=(overlap<=abases.length ? abases.length-overlap : 0);
		//			System.err.print("o");

		for(int i=istart, j=jstart; i<overlap && i<bbases.length && j<abases.length; i++, j++){
			final byte ca=abases[j], cb=bbases[i];
			final byte qa=aqual[j], qb=bqual[i];
			
			if(ca=='N' || cb=='N'){
				//do nothing
			}else{
				assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
				float prob=QualityTools.qualsToProbError(qa, qb);
				if(ca==cb){
					good++;
					goodWeighted+=(1-prob);
				}else{
					bad++;
					badWeighted+=prob;
				}
				expected+=prob;
			}
		}
		if(badWeighted>expected){return 0;}
		return expected;
	}
	
	private static final int countMismatches(Read a, Read b, int insert, int maxMismatches, byte minq){
		final int lengthSum=a.length()+b.length();
		if(insert>=lengthSum){return 0;}
		final int overlap=Tools.min(insert, lengthSum-insert);
		
		int mismatches=0;
		
		
		int start1=(insert>a.length() ? a.length()-overlap : 0);
		int start2=(insert>=b.length() ? 0 : b.length()-overlap);
		
		while(start1<0 || start2<0){start1++; start2++;}
		for(int i=start1, j=start2; i<a.length() && j<b.length(); i++, j++){
			final byte ca=a.bases[i], cb=b.bases[j];
			if(ca!=cb){
				final byte qa=a.quality[i], qb=b.quality[j];
				if(ca=='N' || cb=='N' || qa<minq || qb<minq){
					//do nothing
				}else{
					mismatches++;
					if(mismatches>maxMismatches){break;}
				}
			}
		}
		return mismatches;
	}
	
	protected static int minCoverage(final Read r, final KCountArray kca, final int k, boolean makeCanonical, int cutoff){
		final byte[] bases=r.bases;
		if(bases==null || bases.length<k){return cutoff;}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		
		int len=0;
		long kmer=0;
		int min=cutoff;
		
		for(int i=0; i<bases.length && min>=cutoff; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;

				if(len>=k){
					int cov=kca.read(kmer, k, makeCanonical);
					min=Tools.min(min, cov);
				}
			}
		}
		return min;
	}
	
	protected static int calcMinOverlapByEntropy(byte[] bases, int k, short[] counts, int minscore){
		final int bits=2*k;
		final int mask=~((-1)<<(bits));
		int kmer=0, len=0, ones=0, twos=0;
		
		Arrays.fill(counts, (short)0);
		
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			if(i<bases.length){
				final byte b=bases[j];
				if(!AminoAcid.isFullyDefined(b)){
					len=0;
					kmer=0;
				}else{
					len++;
					final int n=Dedupe.baseToNumber[b];
					kmer=((kmer<<2)|n)&mask;

					if(len>=k){
						counts[kmer]++;
						if(counts[kmer]==1){ones++;}
						else if(counts[kmer]==2){twos++;}
						if(ones*4+twos>=minscore){return i;}
					}
				}
			}
		}
		return bases.length+1;
	}
	
	protected static int MIN_OVERLAP_INSERT=25;
	private static final int BAD_MULT=6;
	private static final int GOOD_MULT_1=8;
	private static final int GOOD_MULT_2=400;

	/** Skip alignment and calculate insert from mapping info */ 
	protected static boolean USE_MAPPING=false;
	protected static final boolean IGNORE_MAPPING_STRAND=false;
	
	protected static boolean verbose=false;
	private static final int FAKE_LEN=4000;
	private static final byte FAKE_Q=20;
	private static final byte[] fakeQual=QualityTools.fakeQuality(FAKE_Q, FAKE_LEN);
	
}
