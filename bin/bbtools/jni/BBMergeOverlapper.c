#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include "jgi_BBMergeOverlapper.h"

// C doesn't have min() or max() so we define our own
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; }) 

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; }) 

#define USE_MAPPING 0
#define MIN_OVERLAP_INSERT 25
#define BAD_MULT 6
#define GOOD_MULT_1 8
#define GOOD_MULT_2 400

jint mateByOverlap(jbyte * a_bases, jint a_bases_length, jbyte * b_bases, jint b_bases_length, jbyte * a_quality, jbyte * b_quality, jint * rvector, jint minOverlap0, const jint minOverlap, jint margin, const jint maxMismatches0, const jint maxMismatches, const jint minq, jint a_insertSizeMapped) {

	minOverlap0=min(max(1, minOverlap0), minOverlap);
	margin=max(margin, 0);
	//if(rvector==NULL){rvector=new jint[5];}
	if(USE_MAPPING){
		rvector[0]=100;
		rvector[1]=20;
		rvector[2]=0;
		rvector[4]=0;
		return a_insertSizeMapped;
	}

	const jbyte *abases=a_bases, *bbases=b_bases;
	const jint alen=a_bases_length, blen=b_bases_length;
	jbyte *aqual=NULL; 
	jbyte *bqual=NULL;
	if(a_quality!=NULL){
		aqual=a_quality;
	}
	if(b_quality!=NULL){
		bqual=b_quality;
	}

	jint bestOverlap=-1;
	jint bestGood=-1;
	jint bestBad=maxMismatches0;

	jboolean ambig=0;
	const jint maxOverlap=alen+blen-max(minOverlap, MIN_OVERLAP_INSERT);

	if(aqual!=NULL && bqual!=NULL) {
		// Precomputed matrix
		const jbyte qualsToPhred[42*42] =
		{0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		 1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		 1, 1, 2, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
		 1, 1, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
		 1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
		 1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,10,10,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,
		 1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9,10,10,11,11,11,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,
		 1, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 9,10,11,11,12,12,12,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
		 1, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9,10,10,11,12,12,13,13,13,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
		 1, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9,10,11,11,12,13,13,14,14,14,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 9,10,11,12,12,13,14,14,14,15,15,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,11,12,13,13,14,14,15,15,16,16,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,10,11,12,13,14,14,15,15,16,16,17,17,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,11,12,13,14,15,15,16,16,17,17,18,18,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,13,14,15,16,17,17,18,18,19,19,20,20,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,14,15,16,17,18,18,19,19,20,20,21,21,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,15,16,17,18,19,19,20,20,21,21,22,22,23,23,23,23,23,23,24,24,24,24,24,24,24,24,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,15,16,17,18,19,20,20,21,21,22,22,23,23,24,24,24,24,24,24,25,25,25,25,25,25,25,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,16,17,18,19,20,21,21,22,22,23,23,24,24,25,25,25,25,25,25,26,26,26,26,26,26,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,17,18,19,20,21,22,22,23,23,24,24,25,25,26,26,26,26,26,26,27,27,27,27,27,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,18,19,20,21,22,23,23,24,24,25,25,26,26,27,27,27,27,27,27,28,28,28,28,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,19,20,21,22,23,24,24,25,25,26,26,27,27,28,28,28,28,28,28,29,29,29,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,20,21,22,23,24,25,25,26,26,27,27,28,28,29,29,29,29,29,29,30,30,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,21,22,23,24,25,26,26,27,27,28,28,29,29,30,30,30,30,30,30,31,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,22,23,24,25,26,27,27,28,28,29,29,30,30,31,31,31,31,31,31,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,23,24,25,26,27,28,28,29,29,30,30,31,31,32,32,32,32,32,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,24,25,26,27,28,29,29,30,30,31,31,32,32,33,33,33,33,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,25,26,27,28,29,30,30,31,31,32,32,33,33,34,34,34,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26,27,28,29,30,31,31,32,32,33,33,34,34,35,35,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,27,28,29,30,31,32,32,33,33,34,34,35,35,36,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,28,29,30,31,32,33,33,34,34,35,35,36,36,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,33,34,34,35,35,36,36,37,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,30,31,32,33,34,35,35,36,36,37,37,
		 1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,31,32,33,34,35,36,36,37,37,38};

		for(jint overlap=max(minOverlap0, 0); overlap<maxOverlap; overlap++){
			jint good=0, bad=0;
			jint istart=(overlap<=alen ? 0 : overlap-alen);
			jint jstart=(overlap<=alen ? alen-overlap : 0);
			{
				const jint iters=min(overlap-istart, min(blen-istart, alen-jstart));
				const jint maxi=istart+iters;
				const jint badlim=bestBad+margin;
				jint iter=0;
				for(jint i1=istart+iter, j1=jstart+iter; i1<maxi && bad<=badlim; i1++, j1++){
					const jbyte ca1=abases[j1], cb1=bbases[i1];
					const jbyte q1=qualsToPhred[aqual[j1]*42+bqual[i1]];

					if(q1<minq || ca1=='N' || cb1=='N'){//do nothing
					}else if(ca1==cb1){good++;}
					else{bad++;}
				}
			}

			if(bad*2<good){
				if(good>minOverlap){
					if(bad<=bestBad){
						if(bad<bestBad || (bad==bestBad && good>bestGood)){
							if(bestBad-bad<margin){ambig=1;}
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
						}else if(bad==bestBad){
							ambig=1;
						}

						if(ambig && bestBad<margin){
							rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
							rvector[1]=bestGood;
							rvector[2]=bestBad;
							rvector[4]=ambig;
							return -1;
						}
					}
				}else if(bad<margin){
					ambig=1;
					rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
					rvector[1]=bestGood;
					rvector[2]=bestBad;
					rvector[4]=ambig;
					return -1;
				}
			}
		}
	}else if (aqual==NULL && bqual==NULL){
		for(jint overlap=max(minOverlap0, 0); overlap<maxOverlap; overlap++){
			jint good=0, bad=0;
			jint istart=(overlap<=alen ? 0 : overlap-alen);
			jint jstart=(overlap<=alen ? alen-overlap : 0);
			{
				const jint iters=min(overlap-istart, min(blen-istart, alen-jstart));
				const jint maxi=istart+iters;
				const jint badlim=bestBad+margin;
				jint iter=0;
				for(jint i1=istart+iter, j1=jstart+iter; i1<maxi && bad<=badlim; i1++, j1++){
					const jbyte ca1=abases[j1], cb1=bbases[i1];
					const jbyte q1=30;

					if(q1<minq || ca1=='N' || cb1=='N'){//do nothing
					}else if(ca1==cb1){good++;}
					else{bad++;}
				}
			}

			if(bad*2<good){
				if(good>minOverlap){
					if(bad<=bestBad){
						if(bad<bestBad || (bad==bestBad && good>bestGood)){
							if(bestBad-bad<margin){ambig=1;}
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
						}else if(bad==bestBad){
							ambig=1;
						}

						if(ambig && bestBad<margin){
							rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
							rvector[1]=bestGood;
							rvector[2]=bestBad;
							rvector[4]=ambig;
							return -1;
						}
					}
				}else if(bad<margin){
					ambig=1;
					rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
					rvector[1]=bestGood;
					rvector[2]=bestBad;
					rvector[4]=ambig;
					return -1;
				}
			}
		}
	}else{
		printf("Not handling the case in the native library where only of the two reads in mateByOverlap() has quality values.\n");
		exit(1);
	}

	if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}

	rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
	rvector[1]=bestGood;
	rvector[2]=bestBad;
	rvector[4]=ambig;

	return (bestOverlap<0 ? -1 : alen+blen-bestOverlap);
}

JNIEXPORT jint JNICALL Java_jgi_BBMergeOverlapper_mateByOverlapJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray a_bases,
							jbyteArray b_bases,
							jbyteArray a_quality,
							jbyteArray b_quality,
							jintArray rvector,
							jint minOverlap0,
							jint minOverlap,
							jint margin,
							jint maxMismatches0,
							jint maxMismatches,
							jint minq,
							jint a_insertSizeMapped
                                                        ) {
   jbyte * ja_quality = NULL;
   jbyte * jb_quality = NULL;

   // Get the size of the read and the reference arrays
   jint a_bases_length = (*env)->GetArrayLength(env, a_bases);
   jint b_bases_length = (*env)->GetArrayLength(env, b_bases);

   // Copy arrays from Java
   jbyte * ja_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_bases, NULL);
   jbyte * jb_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_bases, NULL);
   if(a_quality!=NULL) {ja_quality = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_quality, NULL);}
   if(b_quality!=NULL) {jb_quality = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_quality, NULL);}
   jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

   jint returnVal = mateByOverlap(ja_bases, a_bases_length, jb_bases, b_bases_length, ja_quality, jb_quality, jrvector, minOverlap0, minOverlap, margin, maxMismatches0, maxMismatches, minq, a_insertSizeMapped);

   // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
   (*env)->ReleasePrimitiveArrayCritical(env, a_bases, ja_bases, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, b_bases, jb_bases, JNI_ABORT);
   if(ja_quality!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, a_quality, ja_quality, JNI_ABORT);}
   if(jb_quality!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, b_quality, jb_quality, JNI_ABORT);}
   (*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);

   return returnVal;
}

