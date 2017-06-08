package dna;

import java.io.File;

import jgi.CalcTrueQuality;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.ReadStreamWriter;
import stream.SamLine;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import fileIO.ByteFile;
import fileIO.ReadWrite;

/**
 * @author Brian Bushnell
 * @date Mar 21, 2014
 *
 */
public class Parser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public Parser(){}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean parse(String arg, String a, String b){
		if(isJavaFlag(arg)){return true;}

		if(parseQuality(arg, a, b)){return true;}
		if(parseZip(arg, a, b)){return true;}
		if(parseSam(arg, a, b)){return true;}
		if(parseFasta(arg, a, b)){return true;}
		if(parseCommonStatic(arg, a, b)){return true;}
		if(parseHist(arg, a, b)){return true;}

		if(parseFiles(arg, a, b)){return true;}
		if(parseCommon(arg, a, b)){return true;}
		if(parseTrim(arg, a, b)){return true;}
		if(parseInterleaved(arg, a, b)){return true;}
		if(parseMapping(arg, a, b)){return true;}
		return false;
	}

	public boolean parseCommon(String arg, String a, String b){
		if(a.equals("reads") || a.equals("maxreads")){
			maxReads=Tools.parseKMG(b);
		}else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
			assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
		}else if(a.equals("sampleseed")){
			sampleseed=Long.parseLong(b);
		}else if(a.equals("t") || a.equals("threads")){
			if(b.equalsIgnoreCase("auto")){Shared.setThreads(-1);}
			else{Shared.setThreads(Integer.parseInt(b));}
			System.err.println("Set threads to "+Shared.threads());
		}else if(a.equals("append") || a.equals("app")){
			append=ReadStats.append=Tools.parseBoolean(b);
		}else if(a.equals("overwrite") || a.equals("ow")){
			overwrite=Tools.parseBoolean(b);
		}else if(a.equals("testsize")){
			testsize=Tools.parseBoolean(b);
		}else if(a.equals("breaklen") || a.equals("breaklength")){
			breakLength=Integer.parseInt(b);
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseInterleaved(String arg, String a, String b){
		if(a.equals("testinterleaved")){
			FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
			System.err.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
			setInterleaved=true;
		}else if(a.equals("forceinterleaved")){
			FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
			System.err.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			setInterleaved=true;
		}else if(a.equals("interleaved") || a.equals("int")){
			if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
			else{
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				setInterleaved=true;
			}
		}else if(a.equals("overrideinterleaved")){
			boolean x=Tools.parseBoolean(b);
			ReadStreamByteWriter.ignorePairAssertions=x;
			if(x){setInterleaved=true;}
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseTrim(String arg, String a, String b){
		if(a.equals("forcetrimmod") || a.equals("forcemrimmodulo") || a.equals("ftm")){
			forceTrimModulo=Integer.parseInt(b);
		}else if(a.equals("ftl") || a.equals("forcetrimleft")){
			forceTrimLeft=Integer.parseInt(b);
		}else if(a.equals("ftr") || a.equals("forcetrimright")){
			forceTrimRight=Integer.parseInt(b);
		}else if(a.equals("qtrim") || a.equals("trim")){
			if(b==null || b.length()==0){qtrimRight=qtrimLeft=true;}
			else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrimLeft=true;qtrimRight=false;}
			else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrimLeft=false;qtrimRight=true;}
			else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrimLeft=qtrimRight=true;}
			else if(Character.isDigit(b.charAt(0))){
				trimq=(byte)Integer.parseInt(b);
				qtrimLeft=qtrimRight=true;
			}else{qtrimRight=qtrimLeft=Tools.parseBoolean(b);}
		}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
			if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
				TrimRead.optimalMode=true;
				TrimRead.optimalBias=Float.parseFloat(b);
				assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
			}else{
				TrimRead.optimalMode=Tools.parseBoolean(b);
			}
		}else if(a.equals("trimright") || a.equals("qtrimright")){
			qtrimRight=Tools.parseBoolean(b);
		}else if(a.equals("trimleft") || a.equals("qtrimleft")){
			qtrimLeft=Tools.parseBoolean(b);
		}else if(a.equals("trimq") || a.equals("trimquality")){
			trimq=Byte.parseByte(b);
		}else if(a.equals("trimbadsequence")){
			trimBadSequence=Tools.parseBoolean(b);
		}else if(a.equals("chastityfilter") || a.equals("cf")){
			chastityFilter=Tools.parseBoolean(b);
		}else if(a.equals("requirebothbad") || a.equals("rbb")){
			requireBothBad=Tools.parseBoolean(b);
		}else if(a.equals("removeifeitherbad") || a.equals("rieb")){
			requireBothBad=!Tools.parseBoolean(b);
		}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
			minReadLength=Integer.parseInt(b);
		}else if(a.equals("maxlength") || a.equals("maxlen")){
			maxReadLength=Integer.parseInt(b);
		}else if(a.equals("mingc")){
			minGC=Float.parseFloat(b);
			if(minGC>0){filterGC=true;}
			assert(minGC>=0 && minGC<=1) : "mingc should be a decimal number between 0 and 1, inclusive.";
		}else if(a.equals("maxgc")){
			maxGC=Float.parseFloat(b);
			if(maxGC<1){filterGC=true;}
			assert(minGC>=0 && minGC<=1) : "maxgc should be a decimal number between 0 and 1, inclusive.";
		}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
			minLenFraction=Float.parseFloat(b);
		}else if(a.equals("maxns")){
			maxNs=Integer.parseInt(b);
		}else if(a.equals("minavgquality") || a.equals("maq")){
			minAvgQuality=Byte.parseByte(b);
			//Read.AVERAGE_QUALITY_BY_PROBABILITY=false;
		}else if(a.equals("minavgquality2") || a.equals("maq2")){
			minAvgQuality=Byte.parseByte(b);
			Read.AVERAGE_QUALITY_BY_PROBABILITY=true;
		}else if(a.equals("minavgquality1") || a.equals("maq1")){
			minAvgQuality=Byte.parseByte(b);
			Read.AVERAGE_QUALITY_BY_PROBABILITY=false;
		}else if(a.equals("averagequalitybyprobability") || a.equals("aqbp")){
			Read.AVERAGE_QUALITY_BY_PROBABILITY=Tools.parseBoolean(b);
		}else if(a.equals("mintl") || a.equals("mintrimlen") || a.equals("mintrimlength")){
			minTrimLength=Integer.parseInt(b);
		}else if(a.equals("untrim")){
			untrim=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseFiles(String arg, String a, String b){
		if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
			in1=b;
		}else if(a.equals("in2") || a.equals("input2")){
			in2=b;
		}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
			out1=b;
			setOut=true;
		}else if(a.equals("out2") || a.equals("output2")){
			out2=b;
			setOut=true;
		}else if(a.equals("qfin") || a.equals("qfin1")){
			qfin1=b;
		}else if(a.equals("qfout") || a.equals("qfout1")){
			qfout1=b;
			setOut=true;
		}else if(a.equals("qfin2")){
			qfin2=b;
		}else if(a.equals("qfout2")){
			qfout2=b;
			setOut=true;
		}else if(a.equals("extin")){
			extin=b;
		}else if(a.equals("extout")){
			extout=b;
		}else if(a.equals("outsingle") || a.equals("outs")){
			outsingle=b;
			setOut=true;
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseMapping(String arg, String a, String b){
		if(a.equals("idfilter") || a.equals("identityfilter")){
			idFilter=Float.parseFloat(b);
			if(idFilter>1f){idFilter/=100;}
			assert(idFilter<=1f) : "idfilter should be between 0 and 1.";
		}else if(a.equals("subfilter")){
			subfilter=Integer.parseInt(b);
		}else if(a.equals("delfilter")){
			delfilter=Integer.parseInt(b);
		}else if(a.equals("insfilter")){
			insfilter=Integer.parseInt(b);
		}else if(a.equals("indelfilter")){
			indelfilter=Integer.parseInt(b);
		}else if(a.equals("dellenfilter")){
			dellenfilter=Integer.parseInt(b);
		}else if(a.equals("inslenfilter")){
			inslenfilter=Integer.parseInt(b);
		}else if(a.equals("editfilter")){
			editfilter=Integer.parseInt(b);
		}else if(a.equals("build") || a.equals("genome")){
			build=Integer.parseInt(b);
			Data.GENOME_BUILD=build;
		}else{
			return false;
		}
		return true;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean parseCommonStatic(String arg, String a, String b){
		if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription") || a.equals("trimreaddescriptions")){
			Shared.TRIM_READ_COMMENTS=Tools.parseBoolean(b);
		}else if(a.equals("tuc") || a.equals("touppercase")){
			Read.TO_UPPER_CASE=Tools.parseBoolean(b);
		}else if(a.equals("lctn") || a.equals("lowercaseton")){
			Read.LOWER_CASE_TO_N=Tools.parseBoolean(b);
		}else if(a.equals("changequality") || a.equals("cq")){
			Read.CHANGE_QUALITY=Tools.parseBoolean(b);
		}else if(a.equals("tossbrokenreads") || a.equals("tbr")){
			boolean x=Tools.parseBoolean(b);
			Read.NULLIFY_BROKEN_QUALITY=x;
			ConcurrentReadInputStream.REMOVE_DISCARDED_READS=x;
		}else if(a.equals("bf1")){
			ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
			ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
		}else if(a.equals("bf2")){
			ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
			ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
		}else if(a.equals("usejni") || a.equals("jni")){
			Shared.USE_JNI=Tools.parseBoolean(b);
		}else if(a.equals("usempi") || a.equals("mpi")){
			if(b!=null && Character.isDigit(b.charAt(0))){
				Shared.MPI_NUM_RANKS=Integer.parseInt(b);
				Shared.USE_MPI=Shared.MPI_NUM_RANKS>0;
			}else{
				Shared.USE_MPI=Tools.parseBoolean(b);
			}
		}else if(a.equals("crismpi")){
			Shared.USE_CRISMPI=Tools.parseBoolean(b);
		}else if(a.equals("mpikeepall")){
			Shared.MPI_KEEP_ALL=Tools.parseBoolean(b);
		}else if(a.equals("readbufferlength") || a.equals("readbufferlen")){
			Shared.READ_BUFFER_LENGTH=(int)Tools.parseKMG(b);
		}else if(a.equals("readbufferdata")){
			Shared.READ_BUFFER_MAX_DATA=(int)Tools.parseKMG(b);
		}else if(a.equals("readbuffers")){
			Shared.setBuffers(Integer.parseInt(b));
		}else if(a.equals("rbm") || a.equals("renamebymapping")){
			FASTQ.TAG_CUSTOM=Tools.parseBoolean(b);
		}else if(a.equals("don") || a.equals("deleteoldname")){
			FASTQ.DELETE_OLD_NAME=Tools.parseBoolean(b);
		}else if(a.equals("assertcigar")){
			ReadStreamWriter.ASSERT_CIGAR=Tools.parseBoolean(b);
		}else if(a.equals("verbosesamline")){
			SamLine.verbose=Tools.parseBoolean(b);
		}else if(a.equals("parsecustom") || a.equals("fastqparsecustom")){
			FASTQ.PARSE_CUSTOM=Tools.parseBoolean(b);
			System.err.println("Set FASTQ.PARSE_CUSTOM to "+FASTQ.PARSE_CUSTOM);
		}else{
			return false;
		}
		return true;
	}
	
	public static boolean parseQuality(String arg, String a, String b){
		if(a.equals("ignorebadquality") || a.equals("ibq")){
			FASTQ.IGNORE_BAD_QUALITY=Tools.parseBoolean(b);
		}else if(a.equals("ascii") || a.equals("asciioffset") || a.equals("quality") || a.equals("qual")){
			byte x;
			if(b.equalsIgnoreCase("sanger")){x=33;}
			else if(b.equalsIgnoreCase("illumina")){x=64;}
			else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;}
			else{x=(byte)Integer.parseInt(b);}
			qin=qout=x;
		}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
			byte x;
			if(b.equalsIgnoreCase("sanger")){x=33;}
			else if(b.equalsIgnoreCase("illumina")){x=64;}
			else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=true;}
			else{x=(byte)Integer.parseInt(b);}
			qin=x;
		}else if(a.equals("asciiout") || a.equals("qualityout") || a.equals("qualout") || a.equals("qout")){
			byte x;
			if(b.equalsIgnoreCase("sanger")){x=33;}
			else if(b.equalsIgnoreCase("illumina")){x=64;}
			else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY_OUT=true;}
			else{x=(byte)Integer.parseInt(b);}
			qout=x;
		}else if(a.equals("fakequality") || a.equals("qfake")){
			FASTQ.FAKE_QUAL=Byte.parseByte(b);
		}else if(a.equals("fakefastaqual") || a.equals("fakefastaquality") || a.equals("ffq")){
			if(b==null || b.length()<1){b="f";}
			if(Character.isLetter(b.charAt(0))){
				FastaReadInputStream.FAKE_QUALITY=Tools.parseBoolean(b);
			}else{
				int x=Integer.parseInt(b);
				if(x<1){
					FastaReadInputStream.FAKE_QUALITY=false;
				}else{
					FastaReadInputStream.FAKE_QUALITY=true;
					FastaReadInputStream.FAKE_QUALITY_LEVEL=(byte)Tools.min(x, 50);
				}
			}
		}else if(a.equals("qauto")){
			FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
		}else{
			return false;
		}
		return true;
	}
	
	public static boolean parseHist(String arg, String a, String b){
		if(a.equals("qualityhistogram") || a.equals("qualityhist") || a.equals("qhist")){
			ReadStats.QUAL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=(ReadStats.BQUAL_HIST_FILE!=null || ReadStats.QUAL_HIST_FILE!=null || ReadStats.AVG_QUAL_HIST_FILE!=null || ReadStats.BQUAL_HIST_OVERALL_FILE!=null);
			if(ReadStats.COLLECT_QUALITY_STATS){System.err.println("Set quality histogram output to "+ReadStats.QUAL_HIST_FILE);}
		}else if(a.equals("basequalityhistogram") || a.equals("basequalityhist") || a.equals("bqhist")){
			ReadStats.BQUAL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=(ReadStats.BQUAL_HIST_FILE!=null || ReadStats.QUAL_HIST_FILE!=null || ReadStats.AVG_QUAL_HIST_FILE!=null || ReadStats.BQUAL_HIST_OVERALL_FILE!=null);
			if(ReadStats.BQUAL_HIST_FILE!=null){System.err.println("Set bquality histogram output to "+ReadStats.BQUAL_HIST_FILE);}
		}else if(a.equals("averagequalityhistogram") || a.equals("aqhist")){
			ReadStats.AVG_QUAL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=(ReadStats.BQUAL_HIST_FILE!=null || ReadStats.QUAL_HIST_FILE!=null || ReadStats.AVG_QUAL_HIST_FILE!=null || ReadStats.BQUAL_HIST_OVERALL_FILE!=null);
			if(ReadStats.COLLECT_QUALITY_STATS){System.err.println("Set average quality histogram output to "+ReadStats.AVG_QUAL_HIST_FILE);}
		}else if(a.equals("overallbasequalityhistogram") || a.equals("overallbasequalityhist") || a.equals("obqhist")){
			ReadStats.BQUAL_HIST_OVERALL_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=(ReadStats.BQUAL_HIST_FILE!=null || ReadStats.QUAL_HIST_FILE!=null || ReadStats.AVG_QUAL_HIST_FILE!=null || ReadStats.BQUAL_HIST_OVERALL_FILE!=null);
			if(ReadStats.COLLECT_QUALITY_STATS){System.err.println("Set quality histogram output to "+ReadStats.QUAL_HIST_FILE);}
		}else if(a.equals("matchhistogram") || a.equals("matchhist") || a.equals("mhist")){
			ReadStats.MATCH_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_MATCH_STATS=(ReadStats.MATCH_HIST_FILE!=null);
			if(ReadStats.COLLECT_MATCH_STATS){System.err.println("Set match histogram output to "+ReadStats.MATCH_HIST_FILE);}
		}else if(a.equals("inserthistogram") || a.equals("inserthist") || a.equals("ihist")){
			ReadStats.INSERT_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_INSERT_STATS=(ReadStats.INSERT_HIST_FILE!=null);
			if(ReadStats.COLLECT_INSERT_STATS){System.err.println("Set insert size histogram output to "+ReadStats.INSERT_HIST_FILE);}
		}else if(a.equals("basehistogram") || a.equals("basehist") || a.equals("bhist")){
			ReadStats.BASE_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_BASE_STATS=(ReadStats.BASE_HIST_FILE!=null);
			if(ReadStats.COLLECT_BASE_STATS){System.err.println("Set base content histogram output to "+ReadStats.BASE_HIST_FILE);}
		}else if(a.equals("qualityaccuracyhistogram") || a.equals("qahist")){
			ReadStats.QUAL_ACCURACY_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_ACCURACY=(ReadStats.QUAL_ACCURACY_FILE!=null);
			if(ReadStats.COLLECT_QUALITY_ACCURACY){System.err.println("Set quality accuracy histogram output to "+ReadStats.QUAL_ACCURACY_FILE);}
		}else if(a.equals("indelhistogram") || a.equals("indelhist")){
			ReadStats.INDEL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_INDEL_STATS=(ReadStats.INDEL_HIST_FILE!=null);
			if(ReadStats.COLLECT_INDEL_STATS){System.err.println("Set indel histogram output to "+ReadStats.INDEL_HIST_FILE);}
		}else if(a.equals("errorhistogram") || a.equals("ehist")){
			ReadStats.ERROR_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_ERROR_STATS=(ReadStats.ERROR_HIST_FILE!=null);
			if(ReadStats.COLLECT_ERROR_STATS){System.err.println("Set error histogram output to "+ReadStats.ERROR_HIST_FILE);}
		}else if(a.equals("lengthhistogram") || a.equals("lhist")){
			ReadStats.LENGTH_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_LENGTH_STATS=(ReadStats.LENGTH_HIST_FILE!=null);
			if(ReadStats.COLLECT_LENGTH_STATS){System.err.println("Set length histogram output to "+ReadStats.LENGTH_HIST_FILE);}
		}else if(a.equals("gchistogram") || a.equals("gchist")){
			ReadStats.GC_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_GC_STATS=(ReadStats.GC_HIST_FILE!=null);
			if(ReadStats.COLLECT_GC_STATS){System.err.println("Set GC histogram output to "+ReadStats.GC_HIST_FILE);}
		}else if(a.equals("gcbins") || a.equals("gchistbins")){
			if("auto".equalsIgnoreCase(b)){
				ReadStats.GC_BINS=750;
				ReadStats.GC_BINS_AUTO=true;
			}else{
				ReadStats.GC_BINS=Integer.parseInt(b);
				ReadStats.GC_BINS_AUTO=false;
			}
		}else if(a.equals("identityhistogram") || a.equals("idhist")){
			ReadStats.IDENTITY_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_IDENTITY_STATS=(ReadStats.IDENTITY_HIST_FILE!=null);
			if(ReadStats.COLLECT_IDENTITY_STATS){System.err.println("Set identity histogram output to "+ReadStats.IDENTITY_HIST_FILE);}
		}else if(a.equals("idhistlen") || a.equals("idhistlength") || a.equals("idhistbins") || a.equals("idbins")){
			if("auto".equalsIgnoreCase(b)){
				ReadStats.ID_BINS=750;
				ReadStats.ID_BINS_AUTO=true;
			}else{
				ReadStats.ID_BINS=Integer.parseInt(b);
				ReadStats.ID_BINS_AUTO=false;
			}
		}else{
			return false;
		}
		return true;
	}

	public static boolean parseZip(String arg, String a, String b){
		if(a.equals("ziplevel") || a.equals("zl")){
			int x=Integer.parseInt(b);
			if(x>=0){
				ReadWrite.ZIPLEVEL=Tools.min(x, 9);
			}
		}else if(a.equals("usegzip") || a.equals("gzip")){
			ReadWrite.USE_GZIP=Tools.parseBoolean(b);
		}else if(a.equals("usepigz") || a.equals("pigz")){
			if(b!=null && Character.isDigit(b.charAt(0))){
				int zt=Integer.parseInt(b);
				if(zt<1){ReadWrite.USE_PIGZ=false;}
				else{
					ReadWrite.USE_PIGZ=true;
					ReadWrite.MAX_ZIP_THREADS=zt;
					ReadWrite.ZIP_THREAD_DIVISOR=1;
				}
			}else{ReadWrite.USE_PIGZ=Tools.parseBoolean(b);}
		}else if(a.equals("zipthreaddivisor") || a.equals("ztd")){
			ReadWrite.ZIP_THREAD_DIVISOR=Integer.parseInt(b);
		}else if(a.equals("usegunzip") || a.equals("gunzip")){
			ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
		}else if(a.equals("useunpigz") || a.equals("unpigz")){
			ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}

	public static boolean parseSam(String arg, String a, String b){
		if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
			assert(b!=null) : "The sam flag requires a version number, e.g. 'sam=1.4'";
			SamLine.VERSION=Float.parseFloat(b);
		}else if(a.equals("mdtag") || a.equals("md")){
			SamLine.MAKE_MD_TAG=Tools.parseBoolean(b);
		}else if(a.equals("idtag")){
			SamLine.MAKE_IDENTITY_TAG=Tools.parseBoolean(b);
		}else if(a.equals("xmtag") || a.equals("xm")){
			SamLine.MAKE_XM_TAG=Tools.parseBoolean(b);
		}else if(a.equals("stoptag")){
			SamLine.MAKE_STOP_TAG=Tools.parseBoolean(b);
		}else if(a.equals("scoretag")){
			SamLine.MAKE_SCORE_TAG=Tools.parseBoolean(b);
		}else if(a.equals("sortscaffolds")){
			SamLine.SORT_SCAFFOLDS=Tools.parseBoolean(b);
		}else if(a.equals("customtag")){
			SamLine.MAKE_CUSTOM_TAGS=Tools.parseBoolean(b);
		}else if(a.equals("nhtag")){
			SamLine.MAKE_NH_TAG=Tools.parseBoolean(b);
		}else if(a.equals("keepnames")){
			SamLine.KEEP_NAMES=Tools.parseBoolean(b);
		}else if(a.equals("saa") || a.equals("secondaryalignmentasterisks")){
			SamLine.SECONDARY_ALIGNMENT_ASTERISKS=Tools.parseBoolean(b);
		}else if(a.equals("inserttag")){
			SamLine.MAKE_INSERT_TAG=Tools.parseBoolean(b);
		}else if(a.equals("correctnesstag")){
			SamLine.MAKE_CORRECTNESS_TAG=Tools.parseBoolean(b);
		}else if(a.equals("intronlen") || a.equals("intronlength")){
			SamLine.INTRON_LIMIT=Integer.parseInt(b);
			SamLine.setintron=true;
		}else if(a.equals("suppressheader") || a.equals("noheader")){
			ReadStreamWriter.NO_HEADER=Tools.parseBoolean(b);
		}else if(a.equals("noheadersequences") || a.equals("nhs") || a.equals("suppressheadersequences")){
			ReadStreamWriter.NO_HEADER_SEQUENCES=Tools.parseBoolean(b);
		}else if(a.equals("tophat")){
			if(Tools.parseBoolean(b)){
				SamLine.MAKE_TOPHAT_TAGS=true;
				FastaReadInputStream.FAKE_QUALITY=true;
				FastaReadInputStream.FAKE_QUALITY_LEVEL=40;
				SamLine.MAKE_MD_TAG=true;
			}
		}else if(a.equals("xstag") || a.equals("xs")){
			SamLine.MAKE_XS_TAG=true;
			if(b!=null){
				b=b.toLowerCase();
				if(b.startsWith("fr-")){b=b.substring(3);}
				if(b.equals("ss") || b.equals("secondstrand")){
					SamLine.XS_SECONDSTRAND=true;
				}else if(b.equals("fs") || b.equals("firststrand")){
					SamLine.XS_SECONDSTRAND=false;
				}else if(b.equals("us") || b.equals("unstranded")){
					SamLine.XS_SECONDSTRAND=false;
				}else{
					SamLine.MAKE_XS_TAG=Tools.parseBoolean(b);
				}
			}
			SamLine.setxs=true;
		}else if(parseReadgroup(arg, a, b)){
			//do nothing
		}else{
			return false;
		}
		return true;
	}

	public static boolean parseFasta(String arg, String a, String b){
		if(a.equals("fastareadlen") || a.equals("fastareadlength")){
			FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
			FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
		}else if(a.equals("fastaminread") || a.equals("fastaminlen") || a.equals("fastaminlength")){
			FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
		}else if(a.equals("forcesectionname")){
			FastaReadInputStream.FORCE_SECTION_NAME=Tools.parseBoolean(b);
		}else if(a.equals("fastawrap")){
			Shared.FASTA_WRAP=Integer.parseInt(b);
		}else{
			return false;
		}
		return true;
	}
	
	public static boolean parseQualityAdjust(String arg, String a, String b){
		if(a.equals("q102matrix") || a.equals("q102m")){
			CalcTrueQuality.q102matrix=b;
		}else if(a.equals("qbpmatrix") || a.equals("bqpm")){
			CalcTrueQuality.qbpmatrix=b;
		}else if(a.equals("loadq102")){
			CalcTrueQuality.q102=Tools.parseBoolean(b);
		}else if(a.equals("loadqbp")){
			CalcTrueQuality.qbp=Tools.parseBoolean(b);
		}else if(a.equals("loadq10")){
			CalcTrueQuality.q10=Tools.parseBoolean(b);
		}else if(a.equals("loadq12")){
			CalcTrueQuality.q12=Tools.parseBoolean(b);
		}else if(a.equals("loadqb012")){
			CalcTrueQuality.qb012=Tools.parseBoolean(b);
		}else if(a.equals("loadqb234")){
			CalcTrueQuality.qb234=Tools.parseBoolean(b);
		}else if(a.equals("loadqp")){
			CalcTrueQuality.qp=Tools.parseBoolean(b);
		}else if(a.equals("adjustquality") || a.equals("adjq")){
			TrimRead.ADJUST_QUALITY=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}

	public static boolean isJavaFlag(String arg){
		if(arg==null){return false;}
		if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.startsWith("-Xmn") || arg.equals("-ea") || arg.equals("-da")){return true;}
		if(arg.startsWith("Xmx") || arg.startsWith("Xms") || arg.startsWith("Xmn")){
			return arg.length()>3 && Character.isDigit(arg.charAt(3));
		}
		return false;
	}
	

	/** Return true if the user seems confused */
	public static boolean parseHelp(String[] args){
		if(args==null || args.length==0 || (args.length==1 && args[0]==null)){return true;}
		if(args.length>1){return false;}
		final String s=args[0].toLowerCase();
		return s.equals("-h") || s.equals("-help") || s.equals("--help") 
				|| s.equals("-version") || s.equals("--version") || s.equals("?") || s.equals("-?") || (s.equals("help") && !new File(s).exists());
	}
	
	/** Set SamLine Readgroup Strings */
	public static boolean parseReadgroup(String arg, String a, String b){
		if(a.equals("readgroup") || a.equals("readgroupid") || a.equals("rgid")){
			SamLine.READGROUP_ID=b;
			if(b!=null){SamLine.READGROUP_TAG="RG:Z:"+b;}
		}else if(a.equals("readgroupcn") || a.equals("rgcn")){
			SamLine.READGROUP_CN=b;
		}else if(a.equals("readgroupds") || a.equals("rgds")){
			SamLine.READGROUP_DS=b;
		}else if(a.equals("readgroupdt") || a.equals("rgdt")){
			SamLine.READGROUP_DT=b;
		}else if(a.equals("readgroupfo") || a.equals("rgfo")){
			SamLine.READGROUP_FO=b;
		}else if(a.equals("readgroupks") || a.equals("rgks")){
			SamLine.READGROUP_KS=b;
		}else if(a.equals("readgrouplb") || a.equals("rglb")){
			SamLine.READGROUP_LB=b;
		}else if(a.equals("readgrouppg") || a.equals("rgpg")){
			SamLine.READGROUP_PG=b;
		}else if(a.equals("readgrouppi") || a.equals("rgpi")){
			SamLine.READGROUP_PI=b;
		}else if(a.equals("readgrouppl") || a.equals("rgpl")){
			SamLine.READGROUP_PL=b;
		}else if(a.equals("readgrouppu") || a.equals("rgpu")){
			SamLine.READGROUP_PU=b;
		}else if(a.equals("readgroupsm") || a.equals("rgsm")){
			SamLine.READGROUP_SM=b;
		}else{
			return false;
		}
		return true;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int forceTrimModulo=-1;
	public int forceTrimLeft=-1;
	public int forceTrimRight=-1;
	public int build=1;

	public long maxReads=-1;
	public float samplerate=1f;
	public long sampleseed=-1;

	public boolean qtrimLeft=false;
	public boolean qtrimRight=false;
	
	public byte trimq=6;
	public byte minAvgQuality=0;
	public int maxNs=-1;
	public int minReadLength=0;
	public int maxReadLength=-1;
	public int minTrimLength=-1;
	public float minLenFraction=0;
	public float minGC=0;
	public float maxGC=1;
	public boolean filterGC=false;
	public boolean untrim=false;

	public float idFilter=-1;
	public int subfilter=-1;
	public int delfilter=-1;
	public int insfilter=-1;
	public int indelfilter=-1;
	public int dellenfilter=-1;
	public int inslenfilter=-1;
	public int editfilter=-1;
	
	public int breakLength=0;
	/** Toss pair only if both reads are shorter than limit */ 
	public boolean requireBothBad=false;
	public boolean trimBadSequence=false;
	public boolean chastityFilter=false;

	public boolean overwrite=false;
	public boolean append=false;
	public boolean testsize=false;
	
	public boolean setInterleaved=false;
	
	public String in1=null;
	public String in2=null;
	
	public String qfin1=null;
	public String qfin2=null;

	public String out1=null;
	public String out2=null;
	public String outsingle=null;
	public boolean setOut=false;

	public String qfout1=null;
	public String qfout2=null;
	
	public String extin=null;
	public String extout=null;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static byte qin=-1;
	public static byte qout=-1;
	
}
