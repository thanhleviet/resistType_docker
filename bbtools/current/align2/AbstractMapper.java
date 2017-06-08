package align2;

import java.io.File;
import java.io.PrintStream;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;

import jgi.CalcTrueQuality;
import jgi.CoveragePileup;

import stream.ConcurrentLegacyReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentSolidInputStream;
import stream.FASTQ;
import stream.ConcurrentReadOutputStream;
import stream.RandomReadInputStream3;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SequentialReadInputStream;

import dna.ChromosomeArray;
import dna.Data;
import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;

/**
 * Abstract superclass created from BBMap variants.
 * Handles argument parsing, I/O stream initialization and shutdown,
 * thread management, statistics collection and formatting. 
 * @author Brian Bushnell
 * @date Oct 15, 2013
 *
 */
public abstract class AbstractMapper {

	public AbstractMapper(String[] args){
		if(Shared.COMMAND_LINE==null){
			Shared.COMMAND_LINE=(args==null ? null : args.clone());
			Shared.BBMAP_CLASS=this.getClass().getName();
			int x=Shared.BBMAP_CLASS.lastIndexOf('.');
			if(x>=0){Shared.BBMAP_CLASS=Shared.BBMAP_CLASS.substring(x+1);}
		}
		setDefaults();
		String[] args2=preparse0(args);
		String[] args3=preparse(args2);
		parse(args3);
		postparse(args3);
		setup();
		checkFiles();
	}
	
	void printOptions(){
		sysout.println("For help, please consult readme.txt or run the shellscript with no parameters.");
	}
	
	final void abort(AbstractMapThread[] mtts, String message){
		closeStreams(cris, rosA, rosM, rosU, rosB);
		if(mtts!=null){int x=shutDownThreads(mtts, true);}
		if(message==null){throw new RuntimeException();}
		throw new RuntimeException(message);
	}
	
	/** In megabytes */
	final void adjustThreadsforMemory(long threadMem){
		Runtime rt=Runtime.getRuntime();
		long mmemory=rt.maxMemory()/1000000;
		long tmemory=rt.totalMemory()/1000000;
		long fmemory=rt.freeMemory()/1000000;
		long umemory=tmemory-fmemory;
		long amemory=mmemory-umemory-40;
//		System.err.println("mmemory="+mmemory+", tmemory="+tmemory+", fmemory="+fmemory+", umemory="+umemory+", amemory="+amemory);
		int maxThreads=(int)(amemory/threadMem);
		if(Shared.threads()>maxThreads){
			System.err.println("\nMax Memory = "+mmemory+" MB\nAvailable Memory = "+amemory+" MB");
			if(maxThreads<1){abort(null, "\n\nNot enough memory.  Please run on a node with at least "+((long)((umemory+40+threadMem)*1.15))+" MB.\n");}
			System.err.println("Reducing threads from "+Shared.threads()+" to "+maxThreads+" due to low system memory.");
			Shared.setThreads(maxThreads);
		}
	}
	
	abstract void setDefaults();
	
	abstract String[] preparse(String[] args);
	
	abstract void postparse(String[] args);

	abstract void setup();
	
	abstract void loadIndex();
	
	abstract void processAmbig2();
	
	abstract void testSpeed(String[] args);
	
	abstract void setSemiperfectMode();
	
	abstract void setPerfectMode();

	abstract void printSettings(int k);
	
	private final void parse(String[] args){
		
		
		sysout.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		sysout.println("BBMap version "+Shared.BBMAP_VERSION_STRING);
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer();
		t.start();
		
		Read.TO_UPPER_CASE=true;

		boolean setMaxIndel1=false, setMaxIndel2=false;
		Parser parser=new Parser();
		parser.minTrimLength=minTrimLength;
		
		for(int i=0; i<args.length; i++){
			final String arg=(args[i]==null ? "null" : args[i]);
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
//			System.err.println("Processing "+arg);
			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				if(a.equals("ziplevel") || a.equals("zl")){//Handle conflated term
					ziplevel=Integer.parseInt(b);
				}
			}else if(Parser.parseHist(arg, a, b)){
				//do nothing
			}else if(Parser.parseSam(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(parser.parseMapping(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("printtoerr")){
				if(Tools.parseBoolean(b)){
					sysout=System.err; 
					Data.sysout=System.err;
				}
			}else if(a.equals("colorspace") || a.equals("cs")){
				colorspace=Tools.parseBoolean(b);
				sysout.println("Set colorspace to "+colorspace);
			}else if(a.equals("path") || a.equals("root")){
				Data.setPath(b);
			}else if(a.equals("ref") || a.equals("reference") || a.equals("fasta")){
				reference=b;
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("out")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")){
					outFile=null;
				}else{
					outFile=b;
//					outFile=b.replace('#', '1');
//					outFile2=(b.contains("#") ? b.replace('#', '2') : null);
				}
			}else if(a.equals("out1")){
				outFile=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
				if(outFile==null){
					outFile=null;
				}
			}else if(a.equals("out2")){
				outFile2=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("outmapped") || a.equals("outmapped1")){
				outFileM=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outm2") || a.equals("outmapped2")){
				outFileM2=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("outunmapped") || a.equals("outunmapped1")){
				outFileU=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outu2") || a.equals("outunmapped2")){
				outFileU2=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outb") || a.equals("outb1") || a.equals("outblack") || a.equals("outblack1") || a.equals("outblacklist") || a.equals("outblacklist1")){
				outFileB=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outb2") || a.equals("outblack2") || a.equals("outblacklist2")){
				outFileB2=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("blacklist") && !Data.scaffoldPrefixes){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")){blacklist=null;}
				else{
					if(blacklist==null){blacklist=new ArrayList<String>();}
					if(b.indexOf(',')<0 || new File(b).exists()){blacklist.add(b);}
					else{
						String[] temp=b.split(",");
						for(String tmp : temp){blacklist.add(tmp);}
					}
				}
			}else if(a.startsWith("out_") && b!=null){
				//ignore, it will be processed later
				if(splitterOutputs==null){splitterOutputs=new ArrayList<String>();}
				splitterOutputs.add(b);
			}else if(a.equals("bamscript") || a.equals("bs")){
				bamscript=b;
			}else if(a.equals("local")){
				LOCAL_ALIGN=Tools.parseBoolean(b);
			}else if(a.equals("skipreads")){
				AbstractMapThread.SKIP_INITIAL=Tools.parseKMG(b);
			}else if(a.equals("readlen") || a.equals("length") || a.equals("len")){
				synthReadlen=Integer.parseInt(b);
			}else if(a.equals("kfilter")){
				KFILTER=Integer.parseInt(b);
			}else if(a.equals("msa")){
				MSA_TYPE=b;
			}else if(a.equals("bandwidth") || a.equals("bw")){
				int x=Tools.max(0, Integer.parseInt(b));
				MSA.bandwidth=x;
			}else if(a.equals("bandwidthratio") || a.equals("bwr")){
				float x=Tools.max(0, Float.parseFloat(b));
				MSA.bandwidthRatio=x;
				assert(x>=0) : "Bandwidth ratio should be at least 0.";
			}else if(a.equals("eono") || a.equals("erroronnooutput")){
				ERROR_ON_NO_OUTPUT=Tools.parseBoolean(b);
			}else if(a.equals("log")){
				RefToIndex.LOG=Tools.parseBoolean(b);
			}else if(a.equals("sitesonly") || a.equals("outputsitesonly")){
				outputSitesOnly=Tools.parseBoolean(b);
				sysout.println("Set outputSitesOnly to "+outputSitesOnly);
			}else if(a.equals("discardambiguous") || a.equals("tossambiguous")){
				REMOVE_DUPLICATE_BEST_ALIGNMENTS=Tools.parseBoolean(b);
				sysout.println("Set REMOVE_DUPLICATE_BEST_ALIGNMENTS to "+REMOVE_DUPLICATE_BEST_ALIGNMENTS);
			}else if(a.equals("ambiguous") || a.equals("ambig")){
				if(b==null){
					throw new RuntimeException(arg);
				}else if(b.equalsIgnoreCase("keep") || b.equalsIgnoreCase("best") || b.equalsIgnoreCase("first")){
					ambigMode=AMBIG_BEST;
				}else if(b.equalsIgnoreCase("all")){
					ambigMode=AMBIG_ALL;
				}else if(b.equalsIgnoreCase("random")){
					ambigMode=AMBIG_RANDOM;
				}else if(b.equalsIgnoreCase("toss") || b.equalsIgnoreCase("discard") || b.equalsIgnoreCase("remove")){
					ambigMode=AMBIG_TOSS;
				}else{
					throw new RuntimeException(arg);
				}
//				sysout.println("Set REMOVE_DUPLICATE_BEST_ALIGNMENTS to "+REMOVE_DUPLICATE_BEST_ALIGNMENTS);
			}else if(a.equals("penalizeambiguous") || a.equals("penalizeambig") || a.equals("pambig")){
				AbstractMapThread.PENALIZE_AMBIG=SamLine.PENALIZE_AMBIG=Tools.parseBoolean(b);
			}else if(a.equals("maxsites")){
				MAX_SITESCORES_TO_PRINT=Integer.parseInt(b);
				AbstractMapThread.MAX_TRIM_SITES_TO_RETAIN=Tools.max(MAX_SITESCORES_TO_PRINT*2, AbstractMapThread.MAX_TRIM_SITES_TO_RETAIN);
			}else if(a.equals("maxsites2")){
				AbstractMapThread.MAX_TRIM_SITES_TO_RETAIN=Integer.parseInt(b);
			}else if(a.equals("secondary")){
				PRINT_SECONDARY_ALIGNMENTS=Tools.parseBoolean(b);
				ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS=PRINT_SECONDARY_ALIGNMENTS;
			}else if(a.equals("sssr") || a.equals("secondarysitescoreratio")){
				AbstractMapThread.SECONDARY_SITE_SCORE_RATIO=Float.parseFloat(b);
			}else if(a.equals("ssao") || a.equals("secondarysiteasambiguousonly")){
				AbstractMapThread.PRINT_SECONDARY_ALIGNMENTS_ONLY_FOR_AMBIGUOUS_READS=Tools.parseBoolean(b);
			}else if(a.equals("quickmatch")){
				QUICK_MATCH_STRINGS=Tools.parseBoolean(b);
			}else if(a.equals("ambiguous2") || a.equals("ambig2")){
				if(b==null){
					throw new RuntimeException(arg);
				}else if(b.equalsIgnoreCase("split") || b.equalsIgnoreCase("stream")){
					BBSplitter.AMBIGUOUS2_MODE=BBSplitter.AMBIGUOUS2_SPLIT;
				}else if(b.equalsIgnoreCase("keep") || b.equalsIgnoreCase("best") || b.equalsIgnoreCase("first")){
					BBSplitter.AMBIGUOUS2_MODE=BBSplitter.AMBIGUOUS2_FIRST;
				}else if(b.equalsIgnoreCase("toss") || b.equalsIgnoreCase("discard") || b.equalsIgnoreCase("remove")){
					BBSplitter.AMBIGUOUS2_MODE=BBSplitter.AMBIGUOUS2_TOSS;
				}else if(b.equalsIgnoreCase("random")){
					BBSplitter.AMBIGUOUS2_MODE=BBSplitter.AMBIGUOUS2_RANDOM;
				}else if(b.equalsIgnoreCase("all")){
					BBSplitter.AMBIGUOUS2_MODE=BBSplitter.AMBIGUOUS2_ALL;
				}else{
					throw new RuntimeException(arg);
				}
			}else if(a.equals("forbidselfmapping")){
				FORBID_SELF_MAPPING=Tools.parseBoolean(b);
				sysout.println("Set FORBID_SELF_MAPPING to "+FORBID_SELF_MAPPING);
			}else if(a.equals("match") || a.equals("cigar")){
				if(b!=null){b=b.toLowerCase();}else{b="true";}
				if(b.equals("long") || b.equals("normal")){
					MAKE_MATCH_STRING=true;
					Read.COMPRESS_MATCH_BEFORE_WRITING=false;
//					sysout.println("Writing long match strings.");
				}else if(b.equals("short") || b.equals("compressed")){
					MAKE_MATCH_STRING=true;
					Read.COMPRESS_MATCH_BEFORE_WRITING=true;
//					sysout.println("Writing short match strings.");
				}else{
					MAKE_MATCH_STRING=Tools.parseBoolean(b);
				}

				if(MAKE_MATCH_STRING){
					sysout.println("Cigar strings enabled.");
				}else{
					sysout.println("Cigar strings disabled.");
				}
			}else if(a.equals("semiperfectmode")){
				SEMIPERFECTMODE=Tools.parseBoolean(b);
				if(ziplevel==-1){ziplevel=2;}
			}else if(a.equals("perfectmode")){
				PERFECTMODE=Tools.parseBoolean(b);
				if(ziplevel==-1){ziplevel=2;}
			}else if(a.equals("trimlist")){
				TRIM_LIST=Tools.parseBoolean(b);
			}else if(a.equals("pairedrandom")){
				PAIRED_RANDOM_READS=Tools.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				OUTPUT_ORDERED_READS=Tools.parseBoolean(b);
				sysout.println("Set OUTPUT_ORDERED_READS to "+OUTPUT_ORDERED_READS);
			}else if(a.equals("outputunmapped")){
				DONT_OUTPUT_UNMAPPED_READS=!Tools.parseBoolean(b);
				sysout.println("Set DONT_OUTPUT_UNMAPPED_READS to "+DONT_OUTPUT_UNMAPPED_READS);
			}else if(a.equals("outputblacklisted")){
				DONT_OUTPUT_BLACKLISTED_READS=!Tools.parseBoolean(b);
				sysout.println("Set DONT_OUTPUT_BLACKLISTED_READS to "+DONT_OUTPUT_BLACKLISTED_READS);
			}else if(a.equals("indexloaded")){
				INDEX_LOADED=Tools.parseBoolean(b);
			}else if(a.equals("build") || a.equals("genome") || a.equals("index")){
				build=Integer.parseInt(b);
			}else if(a.equals("minchrom")){
				minChrom=Integer.parseInt(b);
				maxChrom=Tools.max(minChrom, maxChrom);
			}else if(a.equals("maxchrom")){
				maxChrom=Byte.parseByte(b);
				minChrom=Tools.min(minChrom, maxChrom);
			}else if(a.equals("expectedsites")){
				expectedSites=Integer.parseInt(b);
			}else if(a.equals("targetsize")){
				targetGenomeSize=Tools.parseKMG(b);
			}else if(a.equals("fgte")){
				fractionGenomeToExclude=Float.parseFloat(b);
				sysout.println("Set fractionGenomeToExclude to "+String.format("%.4f",fractionGenomeToExclude));
			}else if(a.equals("minratio")){
				MINIMUM_ALIGNMENT_SCORE_RATIO=Float.parseFloat(b);
				sysout.println("Set MINIMUM_ALIGNMENT_SCORE_RATIO to "+String.format("%.3f",MINIMUM_ALIGNMENT_SCORE_RATIO));
				minid=-1;
			}else if(a.equals("minidentity") || a.equals("minid")){
				if(b.lastIndexOf('%')==b.length()-1){minid=Double.parseDouble(b.substring(b.length()-1))/100;}
				else{minid=Double.parseDouble(b);}
				assert(minid>=0 && minid<=100) : "min identity must be between 0 and 1.  Values from 1 to 100 will be assumed percent and divided by 100.";
			}else if(a.equals("rcompmate") || a.equals("reversecomplementmate")){
				rcompMate=Tools.parseBoolean(b);
				sysout.println("Set RCOMP_MATE to "+rcompMate);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				TranslateColorspaceRead.verbose=verbose;
				AbstractIndex.verbose2=verbose;
			}else if(a.equals("verbosestats")){
				if(Character.isDigit(b.charAt(0))){
					verbose_stats=Integer.parseInt(b);
				}else{
					verbose_stats=Tools.parseBoolean(b) ? 9 : 0;
				}
			}else if(a.equals("maxdellen")){
				maxDelLen=Integer.parseInt(b);
			}else if(a.equals("maxinslen")){
				maxInsLen=Integer.parseInt(b);
			}else if(a.equals("maxsublen")){
				maxSubLen=Integer.parseInt(b);
			}else if(a.equals("minqual")){
				minQuality=Byte.parseByte(b);
				midQuality=Tools.max(minQuality, midQuality);
				maxQuality=Tools.max(midQuality, maxQuality);
			}else if(a.equals("midqual")){
				midQuality=Byte.parseByte(b);
				maxQuality=Tools.max(midQuality, maxQuality);
				minQuality=Tools.min(minQuality, midQuality);
			}else if(a.equals("maxqual")){
				maxQuality=Byte.parseByte(b);
				midQuality=Tools.min(maxQuality, midQuality);
				minQuality=Tools.min(minQuality, midQuality);
			}else if(a.equals("matelen") || a.equals("pairlen")){
				int x=Integer.parseInt(b);
				AbstractMapThread.MAX_PAIR_DIST=x;
			}else if(a.equals("s") || a.equals("snps")){
				maxSnps=Integer.parseInt(b);
				baseSnpRate=1;
			}else if(a.equals("u") || a.equals("subs")){
				maxInss=Integer.parseInt(b);
				baseInsRate=1;
			}else if(a.equals("d") || a.equals("dels")){
				maxDels=Integer.parseInt(b);
				baseDelRate=1;
			}else if(a.equals("i") || a.equals("inss")){
				maxSubs=Integer.parseInt(b);
				baseSubRate=1;
			}else if(a.equals("sequentialoverlap")){
				sequentialOverlap=Integer.parseInt(b);
			}else if(a.equals("sequentialstrandalt")){
				sequentialStrandAlt=Tools.parseBoolean(b);
			}else if(a.equals("k") || a.equals("keylen")){
				keylen=Integer.parseInt(b);
			}else if(a.equals("genscaffoldinfo")){
				RefToIndex.genScaffoldInfo=Tools.parseBoolean(b);
			}else if(a.equals("loadscaffolds")){
				Data.LOAD_SCAFFOLDS=Tools.parseBoolean(b);
			}else if(a.equals("autoRefToIndex.chrombits")){
				if("auto".equalsIgnoreCase(b)){RefToIndex.AUTO_CHROMBITS=true;}
				else{RefToIndex.AUTO_CHROMBITS=Tools.parseBoolean(b);}
			}else if(a.equals("RefToIndex.chrombits") || a.equals("cbits")){
				if("auto".equalsIgnoreCase(b)){RefToIndex.AUTO_CHROMBITS=true;}
				else{
					RefToIndex.AUTO_CHROMBITS=false;
					RefToIndex.chrombits=Integer.parseInt(b);
				}
			}else if(a.equals("requirecorrectstrand") || a.equals("rcs")){
				REQUIRE_CORRECT_STRANDS_PAIRS=Tools.parseBoolean(b);
			}else if(a.equals("samestrandpairs") || a.equals("ssp")){
				SAME_STRAND_PAIRS=Tools.parseBoolean(b);
				if(SAME_STRAND_PAIRS){sysout.println("Warning! SAME_STRAND_PAIRS=true mode is not fully tested.");}
			}else if(a.equals("killbadpairs") || a.equals("kbp")){
				KILL_BAD_PAIRS=Tools.parseBoolean(b);
			}else if(a.equals("pairedonly") || a.equals("po")){
				AbstractMapThread.OUTPUT_PAIRED_ONLY=Tools.parseBoolean(b);
			}else if(a.equals("idmodulo") || a.equals("idmod")){
				idmodulo=Integer.parseInt(b);
			}else if(a.equals("minhits") || a.equals("minapproxhits")){
				minApproxHits=Integer.parseInt(b);
			}else if(a.equals("maxindel")){
				maxIndel1=(int)Tools.max(0, Tools.parseKMG(b));
				if(!setMaxIndel2){maxIndel2=2*maxIndel1;}
			}else if(a.equals("maxindel1") || a.equals("maxindelsingle")){
				maxIndel1=(int)Tools.max(0, Tools.parseKMG(b));
				maxIndel2=Tools.max(maxIndel1, maxIndel2);
				setMaxIndel1=true;
			}else if(a.equals("maxindel2") || a.equals("maxindelsum")){
				maxIndel2=(int)Tools.max(0, Tools.parseKMG(b));
				maxIndel1=Tools.min(maxIndel1, maxIndel2);
				setMaxIndel2=true;
			}else if(a.equals("strictmaxindel")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					maxIndel1=(int)Tools.max(0, Tools.parseKMG(b));
					if(!setMaxIndel2){maxIndel2=2*maxIndel1;}
					STRICT_MAX_INDEL=true;
				}else{
					STRICT_MAX_INDEL=Tools.parseBoolean(b);
				}
			}else if(a.equals("padding")){
				SLOW_ALIGN_PADDING=Integer.parseInt(b);
				SLOW_RESCUE_PADDING=SLOW_ALIGN_PADDING;
			}else if(a.equals("rescue")){
				RESCUE=Tools.parseBoolean(b);
			}else if(a.equals("tipsearch")){
				TIP_SEARCH_DIST=Tools.max(0, Integer.parseInt(b));
			}else if(a.equals("dper") || a.equals("dprr")){
				DOUBLE_PRINT_ERROR_RATE=Tools.parseBoolean(b);
			}else if(a.equals("chromgz")){
				Data.CHROMGZ=Tools.parseBoolean(b);
			}else if(a.equals("nodisk")){
				RefToIndex.NODISK=Tools.parseBoolean(b);
			}else if(a.equals("maxchromlen")){
				RefToIndex.maxChromLen=Tools.parseKMG(b);
			}else if(a.equals("minscaf") || a.equals("mincontig")){
				RefToIndex.minScaf=Integer.parseInt(b);
			}else if(a.equals("midpad") || a.equals("interpad")){
				RefToIndex.midPad=Integer.parseInt(b);
			}else if(a.equals("startpad")){
				RefToIndex.startPad=Integer.parseInt(b);
			}else if(a.equals("stoppad")){
				RefToIndex.stopPad=Integer.parseInt(b);
			}else if(a.equals("forceanalyze")){
				forceanalyze=Tools.parseBoolean(b);
			}else if(a.equals("machineoutput") || a.equals("machineout")){
				MACHINE_OUTPUT=Tools.parseBoolean(b);
			}else if(a.equals("showprogress")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					long x=Tools.max(1, Long.parseLong(b));
					ConcurrentReadInputStream.PROGRESS_INCR=x;
					ConcurrentReadInputStream.SHOW_PROGRESS=(x>0);
				}else{
					ConcurrentReadInputStream.SHOW_PROGRESS=Tools.parseBoolean(b);
				}
			}else if(a.equals("scafstats") || a.equals("scaffoldstats")){
				if(b==null && arg.indexOf('=')<0){b="stdout";}
				if(b==null || b.equalsIgnoreCase("false") || b.equalsIgnoreCase("f") || b.equalsIgnoreCase("none") || b.equalsIgnoreCase("null")){
					BBSplitter.TRACK_SCAF_STATS=false;
					BBSplitter.SCAF_STATS_FILE=null;
					sysout.println("No file specified; not tracking scaffold statistics.");
				}else{
					BBSplitter.TRACK_SCAF_STATS=true;
					BBSplitter.SCAF_STATS_FILE=b;
					sysout.println("Scaffold statistics will be written to "+b);
				}
			}else if(a.equals("setstats") || a.equals("refstats")){
				if(b==null && arg.indexOf('=')<0){b="stdout";}
				if(b==null || b.equalsIgnoreCase("false") || b.equalsIgnoreCase("f") || b.equalsIgnoreCase("none") || b.equalsIgnoreCase("null")){
					BBSplitter.TRACK_SET_STATS=false;
					BBSplitter.SET_STATS_FILE=null;
					sysout.println("No file specified; not tracking reference set statistics.");
				}else{
					BBSplitter.TRACK_SET_STATS=true;
					BBSplitter.SET_STATS_FILE=b;
					sysout.println("Reference set statistics will be written to "+b);
				}
			}else if(a.equals("camelwalk")){
				AbstractIndex.USE_CAMELWALK=Tools.parseBoolean(b);
			}else if(a.equals("usequality") || a.equals("uq")){
				AbstractIndex.GENERATE_KEY_SCORES_FROM_QUALITY=AbstractIndex.GENERATE_BASE_SCORES_FROM_QUALITY=Tools.parseBoolean(b);
			}else if(a.equals("ignorequality")){
				AbstractIndex.GENERATE_KEY_SCORES_FROM_QUALITY=AbstractIndex.GENERATE_BASE_SCORES_FROM_QUALITY=!Tools.parseBoolean(b);
			}else if(a.equals("keepbadkeys") || a.equals("kbk")){
				KeyRing.KEEP_BAD_KEYS=Tools.parseBoolean(b);
			}else if(a.equals("usemodulo") || a.equals("um")){
				USE_MODULO=AbstractMapThread.USE_MODULO=IndexMaker4.USE_MODULO=IndexMaker5.USE_MODULO=Tools.parseBoolean(b);
			}else if(a.equals("lowmem") || a.equals("lowram") || a.equals("lowmemory")){
				boolean x=Tools.parseBoolean(b);
				if(x){
					Shared.LOW_MEMORY=true;
					USE_MODULO=AbstractMapThread.USE_MODULO=IndexMaker4.USE_MODULO=IndexMaker5.USE_MODULO=Tools.parseBoolean(b);
				}else{
					Shared.LOW_MEMORY=false;
				}
			}else if(a.equals("coveragestats") || a.equals("covstats")){
				coverageStats=b;
			}else if(a.equals("coverageminscaf") || a.equals("covminscaf")){
				coverageMinScaf=Integer.parseInt(b);
			}else if(a.equals("binnedcoverage") || a.equals("bincov")){
				coverageBinned=b;
			}else if(a.equals("coverage") || a.equals("basecov")){
				coverageBase=b;
			}else if(a.equals("coveragehistogram") || a.equals("covhist")){
				coverageHist=b;
			}else if(a.equals("32bit") || a.equals("32bits") || a.equals("bits32")){
				cov32bit=Tools.parseBoolean(b);
			}else if(a.equals("bitset")){
				covBitset=Tools.parseBoolean(b);
				covSetbs=true;
			}else if(a.equals("arrays")){
				covArrays=Tools.parseBoolean(b);
				covSetbs=true;
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				covNzo=Tools.parseBoolean(b);
			}else if(a.equals("twocolumn")){
				covTwocolumn=Tools.parseBoolean(b);
			}
//			else if(a.equals("uscov") || a.equals("usesecondarycoverage")){
//				covSecondary=Tools.parseBoolean(b);
//			}
			else if(a.equals("ksb") || a.equals("keepshortbins")){
				covKsb=Tools.parseBoolean(b);
			}else if(a.equals("covbinsize")){
				covBinSize=Integer.parseInt(b);
			}else if(a.equals("strandedcoverage") || a.equals("strandedcov") || a.equals("covstranded")){
				covStranded=Tools.parseBoolean(b);
			}else if(a.equals("startcov") || a.equals("covstart")){
				covStartOnly=Tools.parseBoolean(b);
			}else if(a.equals("concisecov")){
				CoveragePileup.CONCISE=Tools.parseBoolean(b);
			}else if(a.equals("printunmappedcount")){
				PRINT_UNMAPPED_COUNT=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+arg);
			}
		}
		
		{//Download parser fields
//			forceTrimLeft=parser.forceTrimLeft;
//			forceTrimRight=parser.forceTrimRight;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			TRIM_QUALITY=parser.trimq;
			AbstractMapThread.MIN_AVERAGE_QUALITY=parser.minAvgQuality;
			AbstractMapThread.MIN_READ_LENGTH=parser.minReadLength;
			AbstractMapThread.MAX_READ_LENGTH=parser.maxReadLength;
//			minLenFraction=parser.minLenFraction;
//			requireBothBad=parser.requireBothBad;
			minTrimLength=parser.minTrimLength;
//			minGC=parser.minGC;
//			maxGC=parser.maxGC;
//			filterGC=parser.filterGC;
//			maxNs=parser.maxNs;
			untrim=parser.untrim;
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=CoveragePileup.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setintron=SamLine.setintron;

			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			IDFILTER=parser.idFilter;
			build=parser.build;
			if(IDFILTER>0){
				if(IDFILTER==1f){PERFECTMODE=true;}
				MAKE_MATCH_STRING=true;
			}
			
			if(parser.subfilter>-1){AbstractMapThread.SUBFILTER=parser.subfilter;}
			if(parser.delfilter>-1){AbstractMapThread.DELFILTER=parser.delfilter;}
			if(parser.insfilter>-1){AbstractMapThread.INSFILTER=parser.insfilter;}
			if(parser.indelfilter>-1){AbstractMapThread.INDELFILTER=parser.indelfilter;}
			if(parser.dellenfilter>-1){AbstractMapThread.DELLENFILTER=parser.dellenfilter;}
			if(parser.inslenfilter>-1){AbstractMapThread.INSLENFILTER=parser.inslenfilter;}
			if(parser.editfilter>-1){AbstractMapThread.EDITFILTER=parser.editfilter;}
		}
		
		gzip=ReadWrite.USE_GZIP;
		gunzip=ReadWrite.USE_GUNZIP;
		pigz=ReadWrite.USE_PIGZ;
		unpigz=ReadWrite.USE_UNPIGZ;
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		ChromosomeArray.UNDEFINED_TO_N=(!INDEX_LOADED);
		
		if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_SPLIT && splitterOutputs!=null){
			ArrayList<String> clone=(ArrayList<String>) splitterOutputs.clone();
			for(String s : clone){
				splitterOutputs.add("AMBIGUOUS_"+s);
			}
		}
	}
	
	private final void checkFiles(){
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){sysout.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(OUTPUT_READS && !Tools.testOutputFiles(overwrite, append, false, outFile, outFile2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outFile+", "+outFile2+"\n");
		}
		
		if(maxReads>0 && maxReads<Long.MAX_VALUE){sysout.println("Max reads: "+maxReads);}
		
		assert(synthReadlen<0 || synthReadlen>=keylen);
	}
	
	private final String[] preparse0(String[] args){
		int nulls=0;
		for(int i=0; i<args.length; i++){
			if(args[i]==null){nulls++;}
			else{
				final String arg=args[i];
				final String[] split=arg.split("=");
				assert(split.length>0) : "\n= symbol must be adjacent to 2 terms, with no spaces.  E.g. 'out=mapped.sam'";
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1].toLowerCase() : null;
				if("null".equalsIgnoreCase(b)){b=null;}
				if(b!=null && (b.equals("stdout") || b.startsWith("stdout."))){
					sysout=System.err;
					Data.sysout=System.err;
				}else if(a.equals("printtoerr")){
					if(Tools.parseBoolean(b)){sysout=System.err; Data.sysout=System.err;}
				}else if(b!=null && (b.equals("stdin") || b.startsWith("stdin."))){
					SYSIN=true;
				}else if(a.equals("fast")){
					fast=Tools.parseBoolean(b);
					if(fast){slow=false;}
					args[i]=null;
					nulls++;
				}else if(a.equals("slow")){
					slow=Tools.parseBoolean(b);
					if(slow){fast=false;}
					args[i]=null;
					nulls++;
				}else if(a.equals("vslow")){
					vslow=Tools.parseBoolean(b);
					if(vslow){fast=false;slow=true;}
					args[i]=null;
					nulls++;
				}
			}
		}
		if(nulls>0){args=Tools.condenseStrict(args);}
		return args;
	}
	
	static final String padPercent(double value, int places){
		String x=String.format("%."+places+"f", value);
		int desired=3+(places<1 ? 0 : 1+places);
		while(x.length()<desired){x=" "+x;}
		return x;
	}
	
	static final String pad(long value, int places){
		String x=""+value;
		while(x.length()<places){x=" "+x;}
		return x;
	}
	
	static final String padPercentMachine(double value, int places){
		String x=String.format("%."+places+"f", value);
		return x;
	}
	

	boolean openStreams(Timer t, String[] args){

		ReadWrite.USE_GZIP=gzip;
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_GUNZIP=gunzip;
		ReadWrite.USE_UNPIGZ=unpigz;
		
		cris=getReadInputStream(in1, in2, qfin1, qfin2);
		final boolean paired=cris.paired();
		cris.setSampleRate(samplerate, sampleseed);
		
		final int buff=(!OUTPUT_ORDERED_READS ? 12 : Tools.max(32, 2*Shared.threads()));
		if(OUTPUT_READS){
			ReadStreamWriter.MINCHROM=minChrom;
			ReadStreamWriter.MAXCHROM=maxChrom;

			if(outFile!=null){
				FileFormat ff1=FileFormat.testOutput(outFile, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				FileFormat ff2=outFile2==null ? null : FileFormat.testOutput(outFile2, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				rosA=ConcurrentReadOutputStream.getStream(ff1, ff2, qfout, qfout2, buff, null, false);
				rosA.start();
				t.stop();
				sysout.println("Started output stream:\t"+t);
				t.start();
			}
			if(outFileM!=null){
				FileFormat ff1=FileFormat.testOutput(outFileM, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				FileFormat ff2=outFileM2==null ? null : FileFormat.testOutput(outFileM2, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				rosM=ConcurrentReadOutputStream.getStream(ff1, ff2, qfoutM, qfoutM2, buff, null, false);
				rosM.start();
				t.stop();
				sysout.println("Started output stream:\t"+t);
				t.start();
			}
			if(outFileU!=null){
				FileFormat ff1=FileFormat.testOutput(outFileU, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				FileFormat ff2=outFileU2==null ? null : FileFormat.testOutput(outFileU2, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				rosU=ConcurrentReadOutputStream.getStream(ff1, ff2, qfoutU, qfoutU2, buff, null, false);
				rosU.start();
				t.stop();
				sysout.println("Started output stream:\t"+t);
				t.start();
			}
			if(outFileB!=null && !Data.scaffoldPrefixes){
				FileFormat ff1=FileFormat.testOutput(outFileB, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				FileFormat ff2=outFileB2==null ? null : FileFormat.testOutput(outFileB2, FileFormat.SAM, 0, 0, true, overwrite, append, OUTPUT_ORDERED_READS);
				rosB=ConcurrentReadOutputStream.getStream(ff1, ff2, qfoutB, qfoutB2, buff, null, false);
				rosB.start();
				t.stop();
				sysout.println("Started output stream:\t"+t);
				t.start();
			}
		}

		if(Data.scaffoldPrefixes){
			BBSplitter.streamTable=BBSplitter.makeOutputStreams(args, OUTPUT_READS, OUTPUT_ORDERED_READS, buff, paired, overwrite, append, false);
			if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_SPLIT){
				BBSplitter.streamTableAmbiguous=BBSplitter.makeOutputStreams(args, OUTPUT_READS, OUTPUT_ORDERED_READS, buff, paired, overwrite, append, true);
			}
		}else{
			BBSplitter.TRACK_SET_STATS=false;
		}

		if(BBSplitter.TRACK_SET_STATS){
			sysout.print("Creating ref-set statistics table: ");
			BBSplitter.makeSetCountTable();
			t.stop();
			sysout.println("   \t"+t);
			t.start();
		}
		if(BBSplitter.TRACK_SCAF_STATS){
			sysout.print("Creating scaffold statistics table:");
			BBSplitter.makeScafCountTable();
			t.stop();
			sysout.println("   \t"+t);
			t.start();
		}

		{
			String syncObj=new String("syncObj");
			synchronized(syncObj){
				System.gc();
				Thread.yield();
//				if(waitForMemoryClear){
					try {syncObj.wait(waitForMemoryClear ? 1000 : 100);} 
					catch (InterruptedException e) {e.printStackTrace();}
//				}
			}

			t.stop();
			sysout.println("Cleared Memory:    \t"+t);
		}
		
		return paired;
	}
	
	static final int shutDownThreads(AbstractMapThread[] mtts, boolean force){
		int broken=0;
		long millis=force ? 500 : 8000;
		for(int i=0; i<mtts.length; i++){
			AbstractMapThread mtt=mtts[i];
			if(mtt==null){broken++;}
			else{
				synchronized(mtt){
					while(mtt.working()){
						State st=mtt.getState();
						if(st==State.TERMINATED){
							if(mtt.working()){
								broken++;
								break;
							}
						}
						try {
							mtt.wait(millis);
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						if(force && mtt.working()){
							mtt.interrupt();
							broken++;
							break;
						}
					}
				}
				if(i==0){
					sysout.print("Detecting finished threads: 0");
				}else{
					sysout.print(", "+i);
				}
			}
		}
		
		if(broken>0){
			System.err.println("\n\n**************************************************************************\n\n" +
					"Warning!  "+broken+" mapping thread"+(broken==1 ? "" : "s")+" did not terminate normally.\n" +
					"Please check the error log; the output may be corrupt or incomplete.\n\n" +
					"**************************************************************************\n\n");
		}
		return broken;
	}
	
	static final boolean closeStreams(ConcurrentReadInputStream cris, ConcurrentReadOutputStream rosA, ConcurrentReadOutputStream rosM, ConcurrentReadOutputStream rosU, ConcurrentReadOutputStream rosB){
		errorState|=ReadWrite.closeStreams(cris, rosA, rosM, rosU, rosB);
		if(BBSplitter.streamTable!=null){
			for(ConcurrentReadOutputStream tros : BBSplitter.streamTable.values()){
				errorState|=ReadWrite.closeStream(tros);
			}
		}
		if(BBSplitter.streamTableAmbiguous!=null){
			for(ConcurrentReadOutputStream tros : BBSplitter.streamTableAmbiguous.values()){
				errorState|=ReadWrite.closeStream(tros);
			}
		}
		return errorState;
	}
	
	static final ConcurrentReadInputStream getReadInputStream(String in1, String in2, String qf1, String qf2){
		
		assert(in1!=null);
		assert(!in1.equalsIgnoreCase(in2)) : in1+", "+in2;
		
		BBIndex.COLORSPACE=colorspace;
		final ConcurrentReadInputStream cris;

		FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, 0, 0, true, true, false);
		FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, 0, 0, true, true, false);
		
		if(ff1.fastq() || ff1.fasta() || ff1.samOrBam() || ff1.scarf() || ff1.bread()){
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, colorspace, false, ff1, ff2, qf1, qf2);
		}else if(ff1.sequential()){
			if(maxReads<0){maxReads=Long.MAX_VALUE;}
//			assert(false) : trials;
			SequentialReadInputStream ris=new SequentialReadInputStream(maxReads, synthReadlen, Tools.max(50, synthReadlen/2), sequentialOverlap, sequentialStrandAlt);
			cris=new ConcurrentLegacyReadInputStream(ris, maxReads);
			
		}else if(ff1.csfasta()){
			colorspace=true;
			BBIndex.COLORSPACE=colorspace;
			
			if(in2!=null){
				cris=new ConcurrentSolidInputStream(in1, qf1, in2, qf2, maxReads);
			}else{
				cris=new ConcurrentSolidInputStream(in1, qf1, maxReads, null);
			}
		}else if(ff1.random()){
			
			useRandomReads=true;
			assert(synthReadlen>0);
			
			RandomReads3.PERFECT_READ_RATIO=PERFECT_READ_RATIO;
			
			RandomReadInputStream3 ris=new RandomReadInputStream3(maxReads, synthReadlen, synthReadlen, 
					maxSnps, maxInss, maxDels, maxSubs,
					baseSnpRate, baseInsRate, baseDelRate, baseSubRate,
					maxInsLen, maxDelLen, maxSubLen,
					minChrom, maxChrom, colorspace, PAIRED_RANDOM_READS,
					minQuality, midQuality, maxQuality);
			cris=new ConcurrentLegacyReadInputStream(ris, maxReads);
		}else{
			throw new RuntimeException("Can't determine read input source: ff1="+ff1+", ff2="+ff2);
		}
		return cris;
	}
	
	
	static void printOutput(final AbstractMapThread[] mtts, final Timer t, final int keylen, final boolean paired, final boolean SKIMMER, final CoveragePileup pile){
		if(MACHINE_OUTPUT){
			printOutput_Machine(mtts, t, keylen, paired, false);
			return;
		}
		
		long readsUsed1=0;
		long readsUsed2=0;
		long lowQualityReadsDiscarded1=0;
		long lowQualityReadsDiscarded2=0;
		long lowQualityBasesDiscarded1=0;
		long lowQualityBasesDiscarded2=0;
		
		long msaIterationsLimited=0;
		long msaIterationsUnlimited=0;

		long basesUsed1=0;
		long basesUsed2=0;
		long keysUsed=0;
		long bothUnmapped=0;
		long bothUnmappedBases=0;
		
		long syntheticReads=0;
		long numMated=0;
		long numMatedBases=0;
		long badPairs=0;
		long badPairBases=0;
		long innerLengthSum=0;
		long outerLengthSum=0;
		long insertSizeSum=0;
		
		long callsToScore=0;
		long callsToExtend=0;
		long initialKeys=0;
		long initialKeyIterations=0;
		long usedKeys=0;
		long usedKeyIterations=0;

		long[] hist_hits=new long[41];
		long[] hist_hits_score=new long[41];
		long[] hist_hits_extend=new long[41];
		
		long initialSiteSum1=0;
		long postTrimSiteSum1=0;
		long postRescueSiteSum1=0;
		long siteSum1=0;
		long topSiteSum1=0;
		
		long matchCountS1=0;
		long matchCountI1=0;
		long matchCountD1=0;
		long matchCountM1=0;
		long matchCountN1=0;
		
		long readCountS1=0;
		long readCountI1=0;
		long readCountD1=0;
		long readCountN1=0;
		long readCountSplice1=0;
		long readCountE1=0;
		
		
		long mapped1=0;
		long mappedRetained1=0;
		long mappedRetainedBases1=0;
		long rescuedP1=0;
		long rescuedM1=0;
		long truePositiveP1=0;
		long truePositiveM1=0;
		long falsePositive1=0;
		long totalCorrectSites1=0;
		long firstSiteCorrectP1=0;
		long firstSiteCorrectM1=0;
		long firstSiteIncorrect1=0;
		long firstSiteCorrectLoose1=0;
		long firstSiteIncorrectLoose1=0;
		long firstSiteCorrectPaired1=0;
		long firstSiteCorrectSolo1=0;
		long firstSiteCorrectRescued1=0;
		long perfectHit1=0; //Highest score is max score
		long uniqueHit1=0; //Only one hit has highest score
		long correctUniqueHit1=0; //unique highest hit on answer site 
		long correctMultiHit1=0;  //non-unique highest hit on answer site (non-skimmer only)
		long correctLowHit1=0;  //hit on answer site, but not highest scorer 
		long noHit1=0;
		long perfectMatch1=0; //Highest slow score is max slow score
		long semiperfectMatch1=0;
		long perfectMatchBases1=0;
		long semiperfectMatchBases1=0;
		long perfectHitCount1=0;
		long semiPerfectHitCount1=0;
		long duplicateBestAlignment1=0;
		long duplicateBestAlignmentBases1=0;
		
		long totalNumCorrect1=0; //Only for skimmer
		long totalNumIncorrect1=0; //Only for skimmer
		long totalNumIncorrectPrior1=0; //Only for skimmer
		long totalNumCapturedAllCorrect1=0; //Only for skimmer
		long totalNumCapturedAllCorrectTop1=0; //Only for skimmer
		long totalNumCapturedAllCorrectOnly1=0; //Only for skimmer

		long initialSiteSum2=0;
		long postTrimSiteSum2=0;
		long postRescueSiteSum2=0;
		long siteSum2=0;
		long topSiteSum2=0;
		
		long mapped2=0;
		long mappedRetained2=0;
		long mappedRetainedBases2=0;
		long rescuedP2=0;
		long rescuedM2=0;
		long truePositiveP2=0;
		long truePositiveM2=0;
		long falsePositive2=0;
		long totalCorrectSites2=0;
		long firstSiteCorrectP2=0;
		long firstSiteCorrectM2=0;
		long firstSiteIncorrect2=0;
		long firstSiteCorrectLoose2=0;
		long firstSiteIncorrectLoose2=0;
		long firstSiteCorrectPaired2=0;
		long firstSiteCorrectSolo2=0;
		long firstSiteCorrectRescued2=0;
		long perfectHit2=0; //Highest score is max score
		long perfectHitCount2=0;
		long semiPerfectHitCount2=0;
		
		long uniqueHit2=0; //Only one hit has highest score
		long correctUniqueHit2=0; //unique highest hit on answer site 
		long correctMultiHit2=0;  //non-unique highest hit on answer site (non-skimmer only)
		long correctLowHit2=0;  //hit on answer site, but not highest scorer 
		long noHit2=0;
		long perfectMatch2=0; //Highest slow score is max slow score
		long semiperfectMatch2=0;
		long perfectMatchBases2=0;
		long semiperfectMatchBases2=0;
		long duplicateBestAlignment2=0;
		long duplicateBestAlignmentBases2=0;
		
		long totalNumCorrect2=0; //Only for skimmer
		long totalNumIncorrect2=0; //Only for skimmer
		long totalNumIncorrectPrior2=0; //Only for skimmer
		long totalNumCapturedAllCorrect2=0; //Only for skimmer
		long totalNumCapturedAllCorrectTop2=0; //Only for skimmer
		long totalNumCapturedAllCorrectOnly2=0; //Only for skimmer
		
		long matchCountS2=0;
		long matchCountI2=0;
		long matchCountD2=0;
		long matchCountM2=0;
		long matchCountN2=0;
		
		long readCountS2=0;
		long readCountI2=0;
		long readCountD2=0;
		long readCountN2=0;
		long readCountSplice2=0;
		long readCountE2=0;
		
		readsUsed1=0;
		for(int i=0; i<mtts.length; i++){
			AbstractMapThread mtt=mtts[i];
			
			if(mtt.msa!=null){
				msaIterationsLimited+=mtt.msa.iterationsLimited;
				msaIterationsUnlimited+=mtt.msa.iterationsUnlimited;
			}
			if(mtt.tcr!=null){
				if(mtt.tcr.msaBS!=null){
					msaIterationsLimited+=mtt.tcr.msaBS.iterationsLimited;
					msaIterationsUnlimited+=mtt.tcr.msaBS.iterationsUnlimited;
				}
				if(mtt.tcr.msaCS!=null){
					msaIterationsLimited+=mtt.tcr.msaCS.iterationsLimited;
					msaIterationsUnlimited+=mtt.tcr.msaCS.iterationsUnlimited;
				}
			}
			
			readsUsed1+=mtt.readsUsed1;
			readsUsed2+=mtt.readsUsed2;
			syntheticReads+=mtt.syntheticReads;
			numMated+=mtt.numMated;
			numMatedBases+=mtt.numMatedBases;
			badPairs+=mtt.badPairs;
			badPairBases+=mtt.badPairBases;
			innerLengthSum+=mtt.innerLengthSum;
			outerLengthSum+=mtt.outerLengthSum;
			insertSizeSum+=mtt.insertSizeSum;
			basesUsed1+=mtt.basesUsed1;
			basesUsed2+=mtt.basesUsed2;
			keysUsed+=mtt.keysUsed;
			bothUnmapped+=mtt.bothUnmapped;
			bothUnmappedBases+=mtt.bothUnmappedBases;
			
			mapped1+=mtt.mapped1;
			mappedRetained1+=mtt.mappedRetained1;
			mappedRetainedBases1+=mtt.mappedRetainedBases1;
			rescuedP1+=mtt.rescuedP1;
			rescuedM1+=mtt.rescuedM1;
			lowQualityReadsDiscarded1+=mtt.lowQualityReadsDiscarded1;
			lowQualityBasesDiscarded1+=mtt.lowQualityBasesDiscarded1;
			truePositiveP1+=mtt.truePositiveP1;
			truePositiveM1+=mtt.truePositiveM1;
			falsePositive1+=mtt.falsePositive1;
//			System.err.println("Adding "+mtt.falsePositive+" false positives -> "+falsePositive);
			totalCorrectSites1+=mtt.totalCorrectSites1;

			firstSiteCorrectP1+=mtt.firstSiteCorrectP1;
			firstSiteCorrectM1+=mtt.firstSiteCorrectM1;
			firstSiteIncorrect1+=mtt.firstSiteIncorrect1;
			firstSiteCorrectLoose1+=mtt.firstSiteCorrectLoose1;
			firstSiteIncorrectLoose1+=mtt.firstSiteIncorrectLoose1;
			firstSiteCorrectPaired1+=mtt.firstSiteCorrectPaired1;
			firstSiteCorrectSolo1+=mtt.firstSiteCorrectSolo1;
			firstSiteCorrectRescued1+=mtt.firstSiteCorrectRescued1;
			
			perfectHit1+=mtt.perfectHit1; //Highest score is max score
			perfectHitCount1+=mtt.perfectHitCount1;
			semiPerfectHitCount1+=mtt.semiPerfectHitCount1;
			uniqueHit1+=mtt.uniqueHit1; //Only one hit has highest score
			correctUniqueHit1+=mtt.correctUniqueHit1; //unique highest hit on answer site 
			correctMultiHit1+=mtt.correctMultiHit1;  //non-unique highest hit on answer site 
			correctLowHit1+=mtt.correctLowHit1;  //hit on answer site, but not highest scorer 
			noHit1+=mtt.noHit1;
			
			totalNumCorrect1+=mtt.totalNumCorrect1; //Skimmer only
			totalNumIncorrect1+=mtt.totalNumIncorrect1; //Skimmer only
			totalNumIncorrectPrior1+=mtt.totalNumIncorrectPrior1; //Skimmer only
			totalNumCapturedAllCorrect1+=mtt.totalNumCapturedAllCorrect1; //Skimmer only
			totalNumCapturedAllCorrectTop1+=mtt.totalNumCapturedAllCorrectTop1; //Skimmer only
			totalNumCapturedAllCorrectOnly1+=mtt.totalNumCapturedAllCorrectOnly1; //Skimmer only
			
			perfectMatch1+=mtt.perfectMatch1; //Highest slow score is max slow score
			semiperfectMatch1+=mtt.semiperfectMatch1; //A semiperfect mapping was found
			perfectMatchBases1+=mtt.perfectMatchBases1;
			semiperfectMatchBases1+=mtt.semiperfectMatchBases1;
			
			duplicateBestAlignment1+=mtt.ambiguousBestAlignment1;
			duplicateBestAlignmentBases1+=mtt.ambiguousBestAlignmentBases1;

			initialSiteSum1+=mtt.initialSiteSum1;
			postTrimSiteSum1+=mtt.postTrimSiteSum1;
			postRescueSiteSum1+=mtt.postRescueSiteSum1;
			siteSum1+=mtt.siteSum1;
			topSiteSum1+=mtt.topSiteSum1;
			
			AbstractIndex index=mtt.index();
			callsToScore+=index.callsToScore;
			callsToExtend+=index.callsToExtendScore;
			initialKeys+=index.initialKeys;
			initialKeyIterations+=index.initialKeyIterations;
			usedKeys+=index.usedKeys;
			usedKeyIterations+=index.usedKeyIterations;
			
			for(int j=0; j<index.hist_hits.length; j++){
				int x=Tools.min(hist_hits.length-1, j);
				hist_hits[x]+=index.hist_hits[j];
				hist_hits_score[x]+=index.hist_hits_score[j];
				hist_hits_extend[x]+=index.hist_hits_extend[j];
			}
			
			matchCountS1+=mtt.matchCountS1;
			matchCountI1+=mtt.matchCountI1;
			matchCountD1+=mtt.matchCountD1;
			matchCountM1+=mtt.matchCountM1;
			matchCountN1+=mtt.matchCountN1;
			
			readCountS1+=mtt.readCountS1;
			readCountI1+=mtt.readCountI1;
			readCountD1+=mtt.readCountD1;
			readCountN1+=mtt.readCountN1;
			readCountSplice1+=mtt.readCountSplice1;
			readCountE1+=mtt.readCountE1;

			mapped2+=mtt.mapped2;
			mappedRetained2+=mtt.mappedRetained2;
			mappedRetainedBases2+=mtt.mappedRetainedBases2;
			rescuedP2+=mtt.rescuedP2;
			rescuedM2+=mtt.rescuedM2;
			lowQualityReadsDiscarded2+=mtt.lowQualityReadsDiscarded2;
			lowQualityBasesDiscarded2+=mtt.lowQualityBasesDiscarded2;
			truePositiveP2+=mtt.truePositiveP2;
			truePositiveM2+=mtt.truePositiveM2;
			falsePositive2+=mtt.falsePositive2;
//			System.err.println("Adding "+mtt.falsePositive+" false positives -> "+falsePositive);
			totalCorrectSites2+=mtt.totalCorrectSites2;

			firstSiteCorrectP2+=mtt.firstSiteCorrectP2;
			firstSiteCorrectM2+=mtt.firstSiteCorrectM2;
			firstSiteIncorrect2+=mtt.firstSiteIncorrect2;
			firstSiteCorrectLoose2+=mtt.firstSiteCorrectLoose2;
			firstSiteIncorrectLoose2+=mtt.firstSiteIncorrectLoose2;
			firstSiteCorrectPaired2+=mtt.firstSiteCorrectPaired2;
			firstSiteCorrectSolo2+=mtt.firstSiteCorrectSolo2;
			firstSiteCorrectRescued2+=mtt.firstSiteCorrectRescued2;
			
			perfectHit2+=mtt.perfectHit2; //Highest score is max score
			perfectHitCount2+=mtt.perfectHitCount2;
			semiPerfectHitCount2+=mtt.semiPerfectHitCount2;
			uniqueHit2+=mtt.uniqueHit2; //Only one hit has highest score
			correctUniqueHit2+=mtt.correctUniqueHit2; //unique highest hit on answer site 
			correctMultiHit2+=mtt.correctMultiHit2;  //non-unique highest hit on answer site 
			correctLowHit2+=mtt.correctLowHit2;  //hit on answer site, but not highest scorer 
			noHit2+=mtt.noHit2;
			
			totalNumCorrect2+=mtt.totalNumCorrect2; //Skimmer only
			totalNumIncorrect2+=mtt.totalNumIncorrect2; //Skimmer only
			totalNumIncorrectPrior2+=mtt.totalNumIncorrectPrior2; //Skimmer only
			totalNumCapturedAllCorrect2+=mtt.totalNumCapturedAllCorrect2; //Skimmer only
			totalNumCapturedAllCorrectTop2+=mtt.totalNumCapturedAllCorrectTop2; //Skimmer only
			totalNumCapturedAllCorrectOnly2+=mtt.totalNumCapturedAllCorrectOnly2; //Skimmer only
			
			perfectMatch2+=mtt.perfectMatch2; //Highest slow score is max slow score
			semiperfectMatch2+=mtt.semiperfectMatch2; //A semiperfect mapping was found
			perfectMatchBases2+=mtt.perfectMatchBases2;
			semiperfectMatchBases2+=mtt.semiperfectMatchBases2;
			
			duplicateBestAlignment2+=mtt.ambiguousBestAlignment2;
			duplicateBestAlignmentBases2+=mtt.ambiguousBestAlignmentBases2;

			initialSiteSum2+=mtt.initialSiteSum2;
			postTrimSiteSum2+=mtt.postTrimSiteSum2;
			postRescueSiteSum2+=mtt.postRescueSiteSum2;
			siteSum2+=mtt.siteSum2;
			topSiteSum2+=mtt.topSiteSum2;
			
			matchCountS2+=mtt.matchCountS2;
			matchCountI2+=mtt.matchCountI2;
			matchCountD2+=mtt.matchCountD2;
			matchCountM2+=mtt.matchCountM2;
			matchCountN2+=mtt.matchCountN2;
			
			readCountS2+=mtt.readCountS2;
			readCountI2+=mtt.readCountI2;
			readCountD2+=mtt.readCountD2;
			readCountN2+=mtt.readCountN2;
			readCountSplice2+=mtt.readCountSplice2;
			readCountE2+=mtt.readCountE2;
			
		}
		maxReads=readsUsed1;
		if(syntheticReads>0){SYNTHETIC=true;}
		
		t.stop();
		long nanos=t.elapsed;
		
		if(verbose_stats>1){
			StringBuilder sb=new StringBuilder(1000);
			sb.append("\n\n###################\n#hits\tcount\tscore\textend\n");
			for(int i=0; i<hist_hits.length; i++){
				sb.append(i+"\t"+hist_hits[i]+"\t"+hist_hits_score[i]+"\t"+hist_hits_extend[i]+"\n");
			}
			try {
				ReadWrite.writeString(sb, "hist_hits.txt", true);
			} catch (Throwable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		final long basesUsed=basesUsed1+basesUsed2;
		
		final double invTrials=1d/maxReads;
		final double invTrials100=100d/maxReads;
		final double invBases100=100d/(basesUsed);
		final double invBases100_1=100d/basesUsed1;
		final double invBases100_2=100d/basesUsed2;
		double invSites100=100d/siteSum1;
		
		final double matedPercent=(numMated*invTrials100);
		final double badPairsPercent=(badPairs*invTrials100);
		final double matedPercentBases=(numMatedBases*invBases100);
		final double badPairsPercentBases=(badPairBases*invBases100);
		final double innerLengthAvg=(innerLengthSum*1d/numMated);
		final double outerLengthAvg=(outerLengthSum*1d/numMated);
		final double insertSizeAvg=(insertSizeSum*1d/numMated);
		
		final double readsPerSecond=((readsUsed1+readsUsed2)*1000000000d)/nanos;
		final double fragsPerSecond=(keysUsed*1000000000d)/nanos;
		final double kiloBasesPerSecond=(basesUsed*1000000d)/nanos;
		
		double perfectHitPercent=(perfectHit1*invTrials100); //Highest score is max score
		double perfectMatchPercent=(perfectMatch1*invTrials100);
		double semiperfectMatchPercent=(semiperfectMatch1*invTrials100);
		double perfectMatchPercentBases=(perfectMatchBases1*invBases100_1);
		double semiperfectMatchPercentBases=(semiperfectMatchBases1*invBases100_1);
		
		double perfectHitCountPercent=perfectHitCount1*invSites100;
		double semiPerfectHitCountPercent=semiPerfectHitCount1*invSites100;
		
		double uniqueHitPercent=(uniqueHit1*invTrials100); //Only one hit has highest score
		double correctUniqueHitPercent=(correctUniqueHit1*invTrials100); //unique highest hit on answer site 
		double correctMultiHitPercent=(correctMultiHit1*invTrials100);  //non-unique highest hit on answer site 
		double correctLowHitPercent=(correctLowHit1*invTrials100);  //hit on answer site, but not highest scorer 
		double ambiguousFound=(duplicateBestAlignment1*invTrials100);
		double ambiguousBasesFound=(duplicateBestAlignmentBases1*invBases100_1);
		double correctHighHitPercent=((correctMultiHit1+correctUniqueHit1)*invTrials100);
		double correctHitPercent=((correctLowHit1+correctMultiHit1+correctUniqueHit1)*invTrials100);

		double mappedB=(mapped1*invTrials100);
		double mappedRetainedB=(mappedRetained1*invTrials100);
		double mappedRetainedBasesB=(mappedRetainedBases1*invBases100_1);
		double rescuedPB=(rescuedP1*invTrials100);
		double rescuedMB=(rescuedM1*invTrials100);
		double falsePositiveB=(firstSiteIncorrect1*invTrials100);
		double falsePositiveLooseB=(firstSiteIncorrectLoose1*invTrials100);
		double truePositivePB=(firstSiteCorrectP1*invTrials100);
		double truePositiveMB=(firstSiteCorrectM1*invTrials100);
		double truePositiveStrict=((firstSiteCorrectP1+firstSiteCorrectM1)*invTrials100);
		double truePositiveLoose=(firstSiteCorrectLoose1*invTrials100);
		double snrStrict=10*Math.log10((firstSiteCorrectM1+firstSiteCorrectP1+0.1)/(firstSiteIncorrect1+0.1));
		double snrLoose=10*Math.log10((firstSiteCorrectLoose1+0.1)/(firstSiteIncorrectLoose1+0.1));
		double truePositivePMRatio=(truePositivePB/truePositiveMB);
		double truePositivePairedB=(firstSiteCorrectPaired1*100d/numMated);
		double truePositiveSoloB=(firstSiteCorrectSolo1*100d/(mappedRetained1-numMated));
		double truePositiveRescuedB=(firstSiteCorrectRescued1*100d/(rescuedP1+rescuedM1));
		double noHitPercent=(noHit1*invTrials100);
		
		long mappedReads, unambiguousReads, mappedBases, unambiguousBases;
		if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
			mappedReads=mappedRetained1+duplicateBestAlignment1;
			unambiguousReads=mappedRetained1;
			mappedBases=mappedRetainedBases1+duplicateBestAlignmentBases1;
			unambiguousBases=mappedRetainedBases1;
		}else{
			mappedReads=mappedRetained1;
			unambiguousReads=mappedRetained1-duplicateBestAlignment1;
			mappedBases=mappedRetainedBases1;
			unambiguousBases=mappedRetainedBases1-duplicateBestAlignmentBases1;
		}
		
		double avgNumCorrect=(SKIMMER ? totalNumCorrect1*invTrials : (totalCorrectSites1/(1d*(truePositiveP1+truePositiveM1))));
		double avgNumIncorrect=totalNumIncorrect1*invTrials; //Skimmer only
		double avgNumIncorrectPrior=totalNumIncorrectPrior1*invTrials; //Skimmer only

		double rateCapturedAllCorrect=totalNumCapturedAllCorrect1*invTrials100; //Skimmer only
		double rateCapturedAllTop=totalNumCapturedAllCorrectTop1*invTrials100; //Skimmer only
		double rateCapturedAllOnly=totalNumCapturedAllCorrectOnly1*invTrials100; //Skimmer only

		double avgCallsToScore=(callsToScore*invTrials);
		double avgCallsToExtendScore=(callsToExtend*invTrials);
		double avgInitialKeys=(initialKeys*1d/initialKeyIterations);
		double avgUsedKeys=(usedKeys*1d/usedKeyIterations);
		
		double avgInitialSites=(initialSiteSum1*invTrials);
		double avgPostTrimSites=(postTrimSiteSum1*invTrials);
		double avgPostRescueSites=(postRescueSiteSum1*invTrials);
		double avgSites=(siteSum1*invTrials);
		double avgPerfectSites=(perfectHitCount1*invTrials);
		double avgSemiPerfectSites=(semiPerfectHitCount1*invTrials);
		double avgTopSites=(topSiteSum1*invTrials);
		double lowQualityReadsDiscardedPercent=(lowQualityReadsDiscarded1*invTrials100);
		double lowQualityBasesDiscardedPercent=(lowQualityBasesDiscarded1*invBases100_1);

		long matchErrors=matchCountS1+matchCountI1+matchCountD1;
		long baseLen=matchCountM1+matchCountI1+matchCountS1+matchCountN1;
		long matchLen=matchCountM1+matchCountI1+matchCountS1+matchCountN1+matchCountD1;
		long refLen=matchCountM1+matchCountS1+matchCountN1+matchCountD1;
		double errorRate=matchErrors*100d/matchLen;
		double matchRate=matchCountM1*100d/matchLen;
		double subRate=matchCountS1*100d/matchLen;
		double delRate=matchCountD1*100d/matchLen;
		double insRate=matchCountI1*100d/matchLen;
		double nRate=matchCountN1*100d/matchLen;
		double readSubRate=readCountS1*100d/mapped1;
		double readDelRate=readCountD1*100d/mapped1;
		double readInsRate=readCountI1*100d/mapped1;
		double readNRate=readCountN1*100d/mapped1;
		double readSpliceRate=readCountSplice1*100d/mapped1;
		double readErrorRate=readCountE1*100d/mapped1;
		
		if(SYNTHETIC && verbose_stats==-1){verbose_stats=Tools.max(verbose_stats,9);}
		
		sysout.println("Reads Used:           \t"+(readsUsed1+readsUsed2)+"\t("+(basesUsed)+" bases)");
		sysout.println();
		
		if(useRandomReads){
			sysout.println("Read Length:          \t"+synthReadlen);
			sysout.println("SNP rate:             \t"+baseSnpRate+"\t(max = "+maxSnps+")");
			sysout.println("INS rate:             \t"+baseInsRate+"\t(max = "+maxInss+", maxLen = "+maxInsLen+")");
			sysout.println("DEL rate:             \t"+baseDelRate+"\t(max = "+maxDels+", maxLen = "+maxDelLen+")");
			sysout.println("SUB rate:             \t"+baseSubRate+"\t(max = "+maxSubs+", maxLen = "+maxSubLen+")");
			sysout.println("minQuality:           \t"+minQuality);
			sysout.println("midQuality:           \t"+midQuality);
			sysout.println("maxQuality:           \t"+maxQuality);
			sysout.println("prefect fraction:     \t"+PERFECT_READ_RATIO);
			sysout.println();
		}

		sysout.println("Mapping:          \t"+t);
		sysout.println(String.format("Reads/sec:       \t%.2f", readsPerSecond));
		sysout.println(String.format("kBases/sec:      \t%.2f", kiloBasesPerSecond));
		double milf=msaIterationsLimited*invTrials;
		double milu=msaIterationsUnlimited*invTrials;
		if(verbose_stats>=1){sysout.println("MSA iterations:   \t"+String.format("%.2fL + %.2fU = %.2f", milf,milu,milf+milu));}
		
		if(paired){
			sysout.println("\n\nPairing data:   \tpct reads\tnum reads \tpct bases\t   num bases");
			sysout.println();
			if(paired){
				sysout.println("mated pairs:     \t"+padPercent(matedPercent,4)+"% \t"+pad(numMated,9)+" \t"+padPercent(matedPercentBases,4)+"% \t"+pad(numMatedBases,12));
				sysout.println("bad pairs:       \t"+padPercent(badPairsPercent,4)+"% \t"+pad(badPairs,9)+" \t"+padPercent(badPairsPercentBases,4)+"% \t"+pad(badPairBases,12));
			}

			sysout.println("insert size avg: \t  "+padPercent(insertSizeAvg,2));
			if(ReadStats.COLLECT_INSERT_STATS){
				if(ReadStats.merged==null){ReadStats.mergeAll();}
				long[] array=ReadStats.merged.insertHist.array;
				double median=Tools.median(array);
				double q1=Tools.percentile(array, 0.25);
				double q3=Tools.percentile(array, 0.75);
				double stdev=Tools.standardDeviationHistogram(array);
				//TODO: Quartiles
				sysout.println("insert 25th %:   \t  "+padPercent(q1,2));
				sysout.println("insert median:   \t  "+padPercent(median,2));
				sysout.println("insert 75th %:   \t  "+padPercent(q3,2));
				sysout.println("insert std dev:  \t  "+padPercent(stdev,2));
				sysout.println("insert mode:     \t  "+Tools.calcMode(array));
			}
			if(verbose_stats>=1){
				sysout.println(String.format("avg inner length:\t  %.2f", innerLengthAvg));
				sysout.println(String.format("avg insert size: \t  %.2f", outerLengthAvg));
			}
		}
		
		if(PRINT_UNMAPPED_COUNT){
			double invReadsUsed100=100.0/(readsUsed1+readsUsed2);
			double invBasesUsed100=100.0/basesUsed;
			double x=bothUnmapped*invReadsUsed100;
			double y=bothUnmappedBases*invBasesUsed100;
			if(!paired){sysout.println();}
			sysout.println("unmapped:        \t"+padPercent(x,4)+"% \t"+pad(bothUnmapped,9)+" \t"+padPercent(y,4)+"% \t"+pad(bothUnmappedBases,12));
		}
		
		sysout.println();
		sysout.println("\nRead 1 data:      \tpct reads\tnum reads \tpct bases\t   num bases");
		if(verbose_stats>=1){
			if(avgInitialKeys>0){sysout.println(String.format("Avg Initial Keys:      \t"+(avgInitialKeys<100?" ":"")+"%.3f", 
					avgInitialKeys));}
			if(avgUsedKeys>0){sysout.println(String.format("Avg Used Keys:         \t"+(avgUsedKeys<100?" ":"")+"%.3f", 
					avgUsedKeys));}
			if(avgCallsToScore>0){sysout.println(String.format("Avg Calls to Score: \t"+(avgCallsToScore<100?" ":"")+"%.3f",
					avgCallsToScore));}
			if(avgCallsToExtendScore>0){sysout.println(String.format("Avg Calls to Extend:\t"+(avgCallsToExtendScore<100?" ":"")+"%.3f", 
					avgCallsToExtendScore));}
			sysout.println();

			sysout.println(String.format("Avg Initial Sites:  \t"+(avgInitialSites<10?" ":"")+"%.3f", avgInitialSites));
			if(TRIM_LIST){sysout.println(String.format("Avg Post-Trim:      \t"+(avgPostTrimSites<10?" ":"")+"%.3f", avgPostTrimSites));}
			if(paired){sysout.println(String.format("Avg Post-Rescue:    \t"+(avgPostRescueSites<10?" ":"")+"%.3f", avgPostRescueSites));}
			sysout.println(String.format("Avg Final Sites:    \t"+(avgSites<10?" ":"")+"%.3f", avgSites));
			sysout.println(String.format("Avg Top Sites:      \t"+(avgTopSites<10?" ":"")+"%.3f", avgTopSites));
			if(verbose_stats>1){
				sysout.println(String.format("Avg Perfect Sites:  \t"+(avgPerfectSites<10?" ":"")+"%.3f    \t"+
						(perfectHitCountPercent<10?" ":"")+"%.3f%%", avgPerfectSites, perfectHitCountPercent));
				sysout.println(String.format("Avg Semiperfect Sites:\t"+(avgSemiPerfectSites<10?" ":"")+"%.3f    \t"+
						(semiPerfectHitCountPercent<10?" ":"")+"%.3f%%", avgSemiPerfectSites, semiPerfectHitCountPercent));
			}

			if(SYNTHETIC){
				sysout.println(String.format("Avg Correct Sites:  \t"+(avgNumCorrect<10?" ":"")+"%.3f", avgNumCorrect));
				if(SKIMMER){	
					sysout.println(String.format("Avg Incorrect Sites:\t"+(avgNumIncorrect<10?" ":"")+"%.3f", avgNumIncorrect));
					sysout.println(String.format("Avg IncorrectP Sites:\t"+(avgNumIncorrectPrior<10?" ":"")+"%.3f", avgNumIncorrectPrior));
				}
			}
		}
		
		sysout.println();
		if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
			double x=ambiguousFound+mappedRetainedB;
			double y=ambiguousBasesFound+mappedRetainedBasesB;
			sysout.println("mapped:          \t"+padPercent(x,4)+"% \t"+pad(mappedReads,9)+" \t"+padPercent(y,4)+"% \t"+pad(mappedBases,12));
			sysout.println("unambiguous:     \t"+padPercent(mappedRetainedB,4)+"% \t"+pad(unambiguousReads,9)+" \t"+padPercent(mappedRetainedBasesB,4)+"% \t"+pad(unambiguousBases,12));
		}else{			
			double x=mappedRetainedB-ambiguousFound;
			double y=mappedRetainedBasesB-ambiguousBasesFound;
			sysout.println("mapped:          \t"+padPercent(mappedRetainedB,4)+"% \t"+pad(mappedReads,9)+" \t"+padPercent(mappedRetainedBasesB,4)+"% \t"+pad(mappedBases,12));
			sysout.println("unambiguous:     \t"+padPercent(x,4)+"% \t"+pad(unambiguousReads,9)+" \t"+padPercent(y,4)+"% \t"+pad(unambiguousBases,12));
		}
		sysout.println("ambiguous:       \t"+padPercent(ambiguousFound,4)+"% \t"+pad(duplicateBestAlignment1,9)+
				" \t"+padPercent(ambiguousBasesFound,4)+"% \t"+pad(duplicateBestAlignmentBases1,12));
		sysout.println("low-Q discards:  \t"+padPercent(lowQualityReadsDiscardedPercent,4)+"% \t"+pad(lowQualityReadsDiscarded1,9)+
				" \t"+padPercent(lowQualityBasesDiscardedPercent,4)+"% \t"+pad(lowQualityBasesDiscarded1,12));
		
		sysout.println();
		sysout.println("perfect best site:\t"+padPercent(perfectMatchPercent,4)+"% \t"+pad(perfectMatch1,9)+
				" \t"+padPercent(perfectMatchPercentBases,4)+"% \t"+pad(perfectMatchBases1,12));
		sysout.println("semiperfect site:\t"+padPercent(semiperfectMatchPercent,4)+"% \t"+pad(semiperfectMatch1,9)+
				" \t"+padPercent(semiperfectMatchPercentBases,4)+"% \t"+pad(semiperfectMatchBases1,12));
		if(paired){
			sysout.println("rescued:         \t"+padPercent(rescuedPB+rescuedMB,4)+"% \t"+pad(rescuedP1+rescuedM1,9));
		}
		
		if(MAKE_MATCH_STRING){
			
			sysout.println();
//			sysout.println("                 \tpct reads\tnum reads \tpct bases\t   num bases");
			sysout.println("Match Rate:      \t      NA \t       NA \t"+padPercent(matchRate,4)+"% \t"+pad(matchCountM1,12));
			sysout.println("Error Rate:      \t"+padPercent(readErrorRate,4)+"% \t"+pad(readCountE1,9)+" \t"+padPercent(errorRate,4)+"% \t"+pad(matchErrors,12));
			sysout.println("Sub Rate:        \t"+padPercent(readSubRate,4)+"% \t"+pad(readCountS1,9)+" \t"+padPercent(subRate,4)+"% \t"+pad(matchCountS1,12));
			sysout.println("Del Rate:        \t"+padPercent(readDelRate,4)+"% \t"+pad(readCountD1,9)+" \t"+padPercent(delRate,4)+"% \t"+pad(matchCountD1,12));
			sysout.println("Ins Rate:        \t"+padPercent(readInsRate,4)+"% \t"+pad(readCountI1,9)+" \t"+padPercent(insRate,4)+"% \t"+pad(matchCountI1,12));
			sysout.println("N Rate:          \t"+padPercent(readNRate,4)+"% \t"+pad(readCountN1,9)+" \t"+padPercent(nRate,4)+"% \t"+pad(matchCountN1,12));
			if(SamLine.INTRON_LIMIT<Integer.MAX_VALUE){
				sysout.println("Splice Rate:     \t"+padPercent(readSpliceRate,4)+"% \t"+pad(readCountSplice1,9)+" \t(splices at least "+SamLine.INTRON_LIMIT+" bp)");
			}
			
			if(DOUBLE_PRINT_ERROR_RATE){
				System.err.println();
				System.err.println(String.format("Match Rate:      \t"+(matchRate<10?" ":"")+"%.4f", matchRate)+"% \t"+matchCountM1);
				System.err.println(String.format("Error Rate:      \t"+(errorRate<10?" ":"")+"%.4f", errorRate)+"% \t"+matchErrors);
				System.err.println(String.format("Sub Rate:        \t"+(subRate<10?" ":"")+"%.4f", subRate)+"% \t"+matchCountS1);
				System.err.println(String.format("Del Rate:        \t"+(delRate<10?" ":"")+"%.4f", delRate)+"% \t"+matchCountD1);
				System.err.println(String.format("Ins Rate:        \t"+(insRate<10?" ":"")+"%.4f", insRate)+"% \t"+matchCountI1);
				System.err.println(String.format("N Rate:          \t"+(nRate<10?" ":"")+"%.4f", nRate)+"% \t"+matchCountN1);
			}
		}
		
		if(SYNTHETIC){
			sysout.println();
			sysout.println("true positive:   \t"+padPercent(truePositiveStrict,4)+"%\t(loose: "+padPercent(truePositiveLoose,4)+"%)");
			sysout.println("false positive:  \t"+padPercent(falsePositiveB,4)+"%\t(loose: "+padPercent(falsePositiveLooseB,4)+"%)");
			sysout.println("false negative:  \t"+padPercent(noHitPercent,4)+"%");
			sysout.println("SNR:             \t"+padPercent(snrStrict,4)+" \t(loose: "+padPercent(snrLoose,4)+")");
			if(verbose_stats>0){
				sysout.println("correctLowHit:   \t"+padPercent(correctLowHitPercent,4)+"%");
				sysout.println(String.format("Plus/Minus ratio:\t %1.4f", truePositivePMRatio));
			}
			
			if(paired){
				sysout.println("correct pairs:   \t"+padPercent(truePositivePairedB,4)+"%\t(of mated)");
				sysout.println("correct singles: \t"+padPercent(truePositiveSoloB,4)+"%");
				sysout.println("correct rescued: \t"+padPercent(truePositiveRescuedB,4)+"%");
			}
			
			if(SKIMMER){
				sysout.println("found all correct:\t"+padPercent(rateCapturedAllCorrect,4)+"%)");
				sysout.println("all correct top:  \t"+padPercent(rateCapturedAllTop,4)+"%)");
				sysout.println("all correct only: \t"+padPercent(rateCapturedAllOnly,4)+"%)");
			}
		}
		
		if(paired){
			
			invSites100=100d/siteSum2;
			
			perfectHitPercent=(perfectHit2*invTrials100); //Highest score is max score
			perfectMatchPercent=(perfectMatch2*invTrials100);
			semiperfectMatchPercent=(semiperfectMatch2*invTrials100);
			perfectMatchPercentBases=(perfectMatchBases2*invBases100_2);
			semiperfectMatchPercentBases=(semiperfectMatchBases2*invBases100_2);
			
			perfectHitCountPercent=perfectHitCount2*invSites100;
			semiPerfectHitCountPercent=semiPerfectHitCount2*invSites100;
			
			uniqueHitPercent=(uniqueHit2*invTrials100); //Only one hit has highest score
			correctUniqueHitPercent=(correctUniqueHit2*invTrials100); //unique highest hit on answer site 
			correctMultiHitPercent=(correctMultiHit2*invTrials100);  //non-unique highest hit on answer site 
			correctLowHitPercent=(correctLowHit2*invTrials100);  //hit on answer site, but not highest scorer 
			ambiguousFound=(duplicateBestAlignment2*invTrials100);
			ambiguousBasesFound=(duplicateBestAlignmentBases2*invBases100_2);
			correctHighHitPercent=((correctMultiHit2+correctUniqueHit2)*invTrials100);
			correctHitPercent=((correctLowHit2+correctMultiHit2+correctUniqueHit2)*invTrials100);

			mappedB=(mapped2*invTrials100);
			mappedRetainedB=(mappedRetained2*invTrials100);
			mappedRetainedBasesB=(mappedRetainedBases2*invBases100_2);
			rescuedPB=(rescuedP2*invTrials100);
			rescuedMB=(rescuedM2*invTrials100);
			falsePositiveB=(firstSiteIncorrect2*invTrials100);
			falsePositiveLooseB=(firstSiteIncorrectLoose2*invTrials100);
			truePositivePB=(firstSiteCorrectP2*invTrials100);
			truePositiveMB=(firstSiteCorrectM2*invTrials100);
			truePositiveStrict=((firstSiteCorrectP2+firstSiteCorrectM2)*invTrials100);
			truePositiveLoose=(firstSiteCorrectLoose2*invTrials100);
			snrStrict=10*Math.log10((firstSiteCorrectM2+firstSiteCorrectP2+0.2)/(firstSiteIncorrect2+0.2));
			snrLoose=10*Math.log10((firstSiteCorrectLoose2+0.2)/(firstSiteIncorrectLoose2+0.2));
			truePositivePMRatio=(truePositivePB/truePositiveMB);
			truePositivePairedB=(firstSiteCorrectPaired2*100d/numMated);
			truePositiveSoloB=(firstSiteCorrectSolo2*100d/(mappedRetained2-numMated));
			truePositiveRescuedB=(firstSiteCorrectRescued2*100d/(rescuedP2+rescuedM2));
			noHitPercent=(noHit2*invTrials100);
			
			if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
				mappedReads=mappedRetained2+duplicateBestAlignment2;
				unambiguousReads=mappedRetained2;
				mappedBases=mappedRetainedBases2+duplicateBestAlignmentBases2;
				unambiguousBases=mappedRetainedBases2;
			}else{
				mappedReads=mappedRetained2;
				unambiguousReads=mappedRetained2-duplicateBestAlignment2;
				mappedBases=mappedRetainedBases2;
				unambiguousBases=mappedRetainedBases2-duplicateBestAlignmentBases2;
			}
			
			avgNumCorrect=(SKIMMER ? totalNumCorrect2*invTrials : (totalCorrectSites2/(2d*(truePositiveP2+truePositiveM2))));
			avgNumIncorrect=totalNumIncorrect2*invTrials; //Skimmer only
			avgNumIncorrectPrior=totalNumIncorrectPrior2*invTrials; //Skimmer only

			rateCapturedAllCorrect=totalNumCapturedAllCorrect2*invTrials100; //Skimmer only
			rateCapturedAllTop=totalNumCapturedAllCorrectTop2*invTrials100; //Skimmer only
			rateCapturedAllOnly=totalNumCapturedAllCorrectOnly2*invTrials100; //Skimmer only

			avgCallsToScore=(callsToScore*invTrials);
			avgCallsToExtendScore=(callsToExtend*invTrials);
			avgInitialKeys=(initialKeys*2d/initialKeyIterations);
			avgUsedKeys=(usedKeys*2d/usedKeyIterations);
			
			avgInitialSites=(initialSiteSum2*invTrials);
			avgPostTrimSites=(postTrimSiteSum2*invTrials);
			avgPostRescueSites=(postRescueSiteSum2*invTrials);
			avgSites=(siteSum2*invTrials);
			avgPerfectSites=(perfectHitCount2*invTrials);
			avgSemiPerfectSites=(semiPerfectHitCount2*invTrials);
			avgTopSites=(topSiteSum2*invTrials);
			lowQualityReadsDiscardedPercent=(lowQualityReadsDiscarded2*invTrials100);
			lowQualityBasesDiscardedPercent=(lowQualityBasesDiscarded2*invBases100_2);

			matchErrors=matchCountS2+matchCountI2+matchCountD2;
			baseLen=matchCountM2+matchCountI2+matchCountS2+matchCountN2;
			matchLen=matchCountM2+matchCountI2+matchCountS2+matchCountN2+matchCountD2;
			refLen=matchCountM2+matchCountS2+matchCountN2+matchCountD2;
			errorRate=matchErrors*100d/matchLen;
			matchRate=matchCountM2*100d/matchLen;
			subRate=matchCountS2*100d/matchLen;
			delRate=matchCountD2*100d/matchLen;
			insRate=matchCountI2*100d/matchLen;
			nRate=matchCountN2*100d/matchLen;
			readSubRate=readCountS2*100d/mapped2;
			readDelRate=readCountD2*100d/mapped2;
			readInsRate=readCountI2*100d/mapped2;
			readNRate=readCountN2*100d/mapped2;
			readSpliceRate=readCountSplice2*100d/mapped2;
			readErrorRate=readCountE2*100d/mapped2;
			
			sysout.println();
			sysout.println("\nRead 2 data:      \tpct reads\tnum reads \tpct bases\t   num bases");
			if(verbose_stats>=1){
				if(avgInitialKeys>0){sysout.println(String.format("Avg Initial Keys:      \t"+(avgInitialKeys<100?" ":"")+"%.3f", 
						avgInitialKeys));}
				if(avgUsedKeys>0){sysout.println(String.format("Avg Used Keys:         \t"+(avgUsedKeys<100?" ":"")+"%.3f", 
						avgUsedKeys));}
				if(avgCallsToScore>0){sysout.println(String.format("Avg Calls to Score: \t"+(avgCallsToScore<100?" ":"")+"%.3f",
						avgCallsToScore));}
				if(avgCallsToExtendScore>0){sysout.println(String.format("Avg Calls to Extend:\t"+(avgCallsToExtendScore<100?" ":"")+"%.3f", 
						avgCallsToExtendScore));}
				sysout.println();

				sysout.println(String.format("Avg Initial Sites:  \t"+(avgInitialSites<10?" ":"")+"%.3f", avgInitialSites));
				if(TRIM_LIST){sysout.println(String.format("Avg Post-Trim:      \t"+(avgPostTrimSites<10?" ":"")+"%.3f", avgPostTrimSites));}
				if(paired){sysout.println(String.format("Avg Post-Rescue:    \t"+(avgPostRescueSites<10?" ":"")+"%.3f", avgPostRescueSites));}
				sysout.println(String.format("Avg Final Sites:    \t"+(avgSites<10?" ":"")+"%.3f", avgSites));
				sysout.println(String.format("Avg Top Sites:      \t"+(avgTopSites<10?" ":"")+"%.3f", avgTopSites));
				if(verbose_stats>1){
					sysout.println(String.format("Avg Perfect Sites:  \t"+(avgPerfectSites<10?" ":"")+"%.3f    \t"+
							(perfectHitCountPercent<10?" ":"")+"%.3f%%", avgPerfectSites, perfectHitCountPercent));
					sysout.println(String.format("Avg Semiperfect Sites:\t"+(avgSemiPerfectSites<10?" ":"")+"%.3f    \t"+
							(semiPerfectHitCountPercent<10?" ":"")+"%.3f%%", avgSemiPerfectSites, semiPerfectHitCountPercent));
				}

				if(SYNTHETIC){
					sysout.println(String.format("Avg Correct Sites:  \t"+(avgNumCorrect<10?" ":"")+"%.3f", avgNumCorrect));
					if(SKIMMER){	
						sysout.println(String.format("Avg Incorrect Sites:\t"+(avgNumIncorrect<10?" ":"")+"%.3f", avgNumIncorrect));
						sysout.println(String.format("Avg IncorrectP Sites:\t"+(avgNumIncorrectPrior<10?" ":"")+"%.3f", avgNumIncorrectPrior));
					}
				}
			}
			
			sysout.println();
			if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
				double x=ambiguousFound+mappedRetainedB;
				double y=ambiguousBasesFound+mappedRetainedBasesB;
				sysout.println("mapped:          \t"+padPercent(x,4)+"% \t"+pad(mappedReads,9)+" \t"+padPercent(y,4)+"% \t"+pad(mappedBases,12));
				sysout.println("unambiguous:     \t"+padPercent(mappedRetainedB,4)+"% \t"+pad(unambiguousReads,9)+" \t"+padPercent(mappedRetainedBasesB,4)+"% \t"+pad(unambiguousBases,12));
			}else{			
				double x=mappedRetainedB-ambiguousFound;
				double y=mappedRetainedBasesB-ambiguousBasesFound;
				sysout.println("mapped:          \t"+padPercent(mappedRetainedB,4)+"% \t"+pad(mappedReads,9)+" \t"+padPercent(mappedRetainedBasesB,4)+"% \t"+pad(mappedBases,12));
				sysout.println("unambiguous:     \t"+padPercent(x,4)+"% \t"+pad(unambiguousReads,9)+" \t"+padPercent(y,4)+"% \t"+pad(unambiguousBases,12));
			}
			sysout.println("ambiguous:       \t"+padPercent(ambiguousFound,4)+"% \t"+pad(duplicateBestAlignment2,9)+
					" \t"+padPercent(ambiguousBasesFound,4)+"% \t"+pad(duplicateBestAlignmentBases2,12));
			sysout.println("low-Q discards:  \t"+padPercent(lowQualityReadsDiscardedPercent,4)+"% \t"+pad(lowQualityReadsDiscarded2,9)+
					" \t"+padPercent(lowQualityBasesDiscardedPercent,4)+"% \t"+pad(lowQualityBasesDiscarded2,12));
			
			sysout.println();
			sysout.println("perfect best site:\t"+padPercent(perfectMatchPercent,4)+"% \t"+pad(perfectMatch2,9)+
					" \t"+padPercent(perfectMatchPercentBases,4)+"% \t"+pad(perfectMatchBases2,12));
			sysout.println("semiperfect site:\t"+padPercent(semiperfectMatchPercent,4)+"% \t"+pad(semiperfectMatch2,9)+
					" \t"+padPercent(semiperfectMatchPercentBases,4)+"% \t"+pad(semiperfectMatchBases2,12));
			if(paired){
				sysout.println("rescued:         \t"+padPercent(rescuedPB+rescuedMB,4)+"% \t"+pad(rescuedP2+rescuedM2,9));
			}
			
			if(MAKE_MATCH_STRING){
				
				sysout.println();
//				sysout.println("                 \tpct reads\tnum reads \tpct bases\t   num bases");
				sysout.println("Match Rate:      \t      NA \t       NA \t"+padPercent(matchRate,4)+"% \t"+pad(matchCountM2,12));
				sysout.println("Error Rate:      \t"+padPercent(readErrorRate,4)+"% \t"+pad(readCountE2,9)+" \t"+padPercent(errorRate,4)+"% \t"+pad(matchErrors,12));
				sysout.println("Sub Rate:        \t"+padPercent(readSubRate,4)+"% \t"+pad(readCountS2,9)+" \t"+padPercent(subRate,4)+"% \t"+pad(matchCountS2,12));
				sysout.println("Del Rate:        \t"+padPercent(readDelRate,4)+"% \t"+pad(readCountD2,9)+" \t"+padPercent(delRate,4)+"% \t"+pad(matchCountD2,12));
				sysout.println("Ins Rate:        \t"+padPercent(readInsRate,4)+"% \t"+pad(readCountI2,9)+" \t"+padPercent(insRate,4)+"% \t"+pad(matchCountI2,12));
				sysout.println("N Rate:          \t"+padPercent(readNRate,4)+"% \t"+pad(readCountN2,9)+" \t"+padPercent(nRate,4)+"% \t"+pad(matchCountN2,12));
				if(SamLine.INTRON_LIMIT<Integer.MAX_VALUE){
					sysout.println("Splice Rate:     \t"+padPercent(readSpliceRate,4)+"% \t"+pad(readCountSplice2,9)+" \t(splices at least "+SamLine.INTRON_LIMIT+" bp)");
				}
				
				if(DOUBLE_PRINT_ERROR_RATE){
					System.err.println();
					System.err.println(String.format("Match Rate:      \t"+(matchRate<10?" ":"")+"%.4f", matchRate)+"% \t"+matchCountM2);
					System.err.println(String.format("Error Rate:      \t"+(errorRate<10?" ":"")+"%.4f", errorRate)+"% \t"+matchErrors);
					System.err.println(String.format("Sub Rate:        \t"+(subRate<10?" ":"")+"%.4f", subRate)+"% \t"+matchCountS2);
					System.err.println(String.format("Del Rate:        \t"+(delRate<10?" ":"")+"%.4f", delRate)+"% \t"+matchCountD2);
					System.err.println(String.format("Ins Rate:        \t"+(insRate<10?" ":"")+"%.4f", insRate)+"% \t"+matchCountI2);
					System.err.println(String.format("N Rate:          \t"+(nRate<10?" ":"")+"%.4f", nRate)+"% \t"+matchCountN2);
				}
			}
			
			if(SYNTHETIC){
				sysout.println();
				sysout.println("true positive:   \t"+padPercent(truePositiveStrict,4)+"%\t(loose: "+padPercent(truePositiveLoose,4)+"%)");
				sysout.println("false positive:  \t"+padPercent(falsePositiveB,4)+"%\t(loose: "+padPercent(falsePositiveLooseB,4)+"%)");
				sysout.println("false negative:  \t"+padPercent(noHitPercent,4)+"%");
				sysout.println("SNR:             \t"+padPercent(snrStrict,4)+" \t(loose: "+padPercent(snrLoose,4)+")");
				if(verbose_stats>0){
					sysout.println("correctLowHit:   \t"+padPercent(correctLowHitPercent,4)+"%");
					sysout.println(String.format("Plus/Minus ratio:\t %2.4f", truePositivePMRatio));
				}
				
				if(paired){
					sysout.println("correct pairs:   \t"+padPercent(truePositivePairedB,4)+"%\t(of mated)");
					sysout.println("correct singles: \t"+padPercent(truePositiveSoloB,4)+"%");
					sysout.println("correct rescued: \t"+padPercent(truePositiveRescuedB,4)+"%");
				}
				
				if(SKIMMER){
					sysout.println("found all correct:\t"+padPercent(rateCapturedAllCorrect,4)+"%)");
					sysout.println("all correct top:  \t"+padPercent(rateCapturedAllTop,4)+"%)");
					sysout.println("all correct only: \t"+padPercent(rateCapturedAllOnly,4)+"%)");
				}
			}
		}
		
		if(BBSplitter.TRACK_SCAF_STATS){
			BBSplitter.printCounts(BBSplitter.SCAF_STATS_FILE, BBSplitter.scafCountTable, true, readsUsed1+readsUsed2);
		}
		
		if(BBSplitter.TRACK_SET_STATS){
			BBSplitter.printCounts(BBSplitter.SET_STATS_FILE, BBSplitter.setCountTable, true, readsUsed1+readsUsed2);
		}
		
		ReadStats.writeAll(paired);
		if(pile!=null){
			CoveragePileup.overwrite=overwrite;
			CoveragePileup.append=append;
			pile.printOutput();
		}
		
		assert(!CALC_STATISTICS || truePositiveP1+truePositiveM1+falsePositive1+noHit1+lowQualityReadsDiscarded1==maxReads) : 
			"\nThe number of reads out does not add up to the number of reads in.\nThis may indicate that a mapping thread crashed.\n"+
			truePositiveP1+"+"+truePositiveM1+"+"+falsePositive1+"+"+noHit1+"+"+lowQualityReadsDiscarded1+" = "+
			(truePositiveP1+truePositiveM1+falsePositive1+noHit1+lowQualityReadsDiscarded1)+" != "+maxReads;
		if(!SKIMMER){
			assert(!CALC_STATISTICS || truePositiveP1+truePositiveM1==correctLowHit1+correctMultiHit1+correctUniqueHit1);
		}else{
			assert(!CALC_STATISTICS || truePositiveP1+truePositiveM1==correctLowHit1+correctUniqueHit1);
		}
	}
	
	
	static void printOutput_Machine(final AbstractMapThread[] mtts, final Timer t, final int keylen, final boolean paired, final boolean SKIMMER){
		
		long readsUsed1=0;
		long readsUsed2=0;
		long lowQualityReadsDiscarded1=0;
		long lowQualityReadsDiscarded2=0;
		long lowQualityBasesDiscarded1=0;
		long lowQualityBasesDiscarded2=0;
		
		long msaIterationsLimited=0;
		long msaIterationsUnlimited=0;

		long basesUsed1=0;
		long basesUsed2=0;
		long basesAtQuickmap=0;
		long keysUsed=0;
		long bothUnmapped=0;
		long bothUnmappedBases=0;
		
		long syntheticReads=0;
		long numMated=0;
		long badPairs=0;
		long innerLengthSum=0;
		long outerLengthSum=0;
		long insertSizeSum=0;
		
		long callsToScore=0;
		long callsToExtend=0;
		long initialKeys=0;
		long initialKeyIterations=0;
		long usedKeys=0;
		long usedKeyIterations=0;

		long[] hist_hits=new long[41];
		long[] hist_hits_score=new long[41];
		long[] hist_hits_extend=new long[41];
		
		long initialSiteSum1=0;
		long postTrimSiteSum1=0;
		long postRescueSiteSum1=0;
		long siteSum1=0;
		long topSiteSum1=0;
		
		long matchCountS1=0;
		long matchCountI1=0;
		long matchCountD1=0;
		long matchCountM1=0;
		long matchCountN1=0;
		
		
		long mapped1=0;
		long mappedRetained1=0;
		long rescuedP1=0;
		long rescuedM1=0;
		long truePositiveP1=0;
		long truePositiveM1=0;
		long falsePositive1=0;
		long totalCorrectSites1=0;
		long firstSiteCorrectP1=0;
		long firstSiteCorrectM1=0;
		long firstSiteIncorrect1=0;
		long firstSiteCorrectLoose1=0;
		long firstSiteIncorrectLoose1=0;
		long firstSiteCorrectPaired1=0;
		long firstSiteCorrectSolo1=0;
		long firstSiteCorrectRescued1=0;
		long perfectHit1=0; //Highest score is max score
		long uniqueHit1=0; //Only one hit has highest score
		long correctUniqueHit1=0; //unique highest hit on answer site 
		long correctMultiHit1=0;  //non-unique highest hit on answer site (non-skimmer only)
		long correctLowHit1=0;  //hit on answer site, but not highest scorer 
		long noHit1=0;
		long perfectMatch1=0; //Highest slow score is max slow score
		long semiperfectMatch1=0;
		long perfectMatchBases1=0;
		long semiperfectMatchBases1=0;
		long perfectHitCount1=0;
		long semiPerfectHitCount1=0;
		long duplicateBestAlignment1=0;
		
		long totalNumCorrect1=0; //Only for skimmer
		long totalNumIncorrect1=0; //Only for skimmer
		long totalNumIncorrectPrior1=0; //Only for skimmer
		long totalNumCapturedAllCorrect1=0; //Only for skimmer
		long totalNumCapturedAllCorrectTop1=0; //Only for skimmer
		long totalNumCapturedAllCorrectOnly1=0; //Only for skimmer

		long initialSiteSum2=0;
		long postTrimSiteSum2=0;
		long postRescueSiteSum2=0;
		long siteSum2=0;
		long topSiteSum2=0;
		
		long mapped2=0;
		long mappedRetained2=0;
		long rescuedP2=0;
		long rescuedM2=0;
		long truePositiveP2=0;
		long truePositiveM2=0;
		long falsePositive2=0;
		long totalCorrectSites2=0;
		long firstSiteCorrectP2=0;
		long firstSiteCorrectM2=0;
		long firstSiteIncorrect2=0;
		long firstSiteCorrectLoose2=0;
		long firstSiteIncorrectLoose2=0;
		long firstSiteCorrectPaired2=0;
		long firstSiteCorrectSolo2=0;
		long firstSiteCorrectRescued2=0;
		long perfectHit2=0; //Highest score is max score
		long perfectHitCount2=0;
		long semiPerfectHitCount2=0;
		
		long uniqueHit2=0; //Only one hit has highest score
		long correctUniqueHit2=0; //unique highest hit on answer site 
		long correctMultiHit2=0;  //non-unique highest hit on answer site (non-skimmer only)
		long correctLowHit2=0;  //hit on answer site, but not highest scorer 
		long noHit2=0;
		long perfectMatch2=0; //Highest slow score is max slow score
		long semiperfectMatch2=0;
		long perfectMatchBases2=0;
		long semiperfectMatchBases2=0;
		long duplicateBestAlignment2=0;
		
		long totalNumCorrect2=0; //Only for skimmer
		long totalNumIncorrect2=0; //Only for skimmer
		long totalNumIncorrectPrior2=0; //Only for skimmer
		long totalNumCapturedAllCorrect2=0; //Only for skimmer
		long totalNumCapturedAllCorrectTop2=0; //Only for skimmer
		long totalNumCapturedAllCorrectOnly2=0; //Only for skimmer
		
		long matchCountS2=0;
		long matchCountI2=0;
		long matchCountD2=0;
		long matchCountM2=0;
		long matchCountN2=0;
		
		readsUsed1=0;
		for(int i=0; i<mtts.length; i++){
			AbstractMapThread mtt=mtts[i];
			
			if(mtt.msa!=null){
				msaIterationsLimited+=mtt.msa.iterationsLimited;
				msaIterationsUnlimited+=mtt.msa.iterationsUnlimited;
			}
			if(mtt.tcr!=null){
				if(mtt.tcr.msaBS!=null){
					msaIterationsLimited+=mtt.tcr.msaBS.iterationsLimited;
					msaIterationsUnlimited+=mtt.tcr.msaBS.iterationsUnlimited;
				}
				if(mtt.tcr.msaCS!=null){
					msaIterationsLimited+=mtt.tcr.msaCS.iterationsLimited;
					msaIterationsUnlimited+=mtt.tcr.msaCS.iterationsUnlimited;
				}
			}
			
			readsUsed1+=mtt.readsUsed1;
			readsUsed2+=mtt.readsUsed2;
			syntheticReads+=mtt.syntheticReads;
			numMated+=mtt.numMated;
			badPairs+=mtt.badPairs;
			innerLengthSum+=mtt.innerLengthSum;
			outerLengthSum+=mtt.outerLengthSum;
			insertSizeSum+=mtt.insertSizeSum;
			basesUsed1+=mtt.basesUsed1;
			basesUsed2+=mtt.basesUsed2;
			keysUsed+=mtt.keysUsed;
			bothUnmapped+=mtt.bothUnmapped;
			bothUnmappedBases+=mtt.bothUnmappedBases;
			
			mapped1+=mtt.mapped1;
			mappedRetained1+=mtt.mappedRetained1;
			rescuedP1+=mtt.rescuedP1;
			rescuedM1+=mtt.rescuedM1;
			lowQualityReadsDiscarded1+=mtt.lowQualityReadsDiscarded1;
			truePositiveP1+=mtt.truePositiveP1;
			truePositiveM1+=mtt.truePositiveM1;
			falsePositive1+=mtt.falsePositive1;
//			System.err.println("Adding "+mtt.falsePositive+" false positives -> "+falsePositive);
			totalCorrectSites1+=mtt.totalCorrectSites1;

			firstSiteCorrectP1+=mtt.firstSiteCorrectP1;
			firstSiteCorrectM1+=mtt.firstSiteCorrectM1;
			firstSiteIncorrect1+=mtt.firstSiteIncorrect1;
			firstSiteCorrectLoose1+=mtt.firstSiteCorrectLoose1;
			firstSiteIncorrectLoose1+=mtt.firstSiteIncorrectLoose1;
			firstSiteCorrectPaired1+=mtt.firstSiteCorrectPaired1;
			firstSiteCorrectSolo1+=mtt.firstSiteCorrectSolo1;
			firstSiteCorrectRescued1+=mtt.firstSiteCorrectRescued1;
			
			perfectHit1+=mtt.perfectHit1; //Highest score is max score
			perfectHitCount1+=mtt.perfectHitCount1;
			semiPerfectHitCount1+=mtt.semiPerfectHitCount1;
			uniqueHit1+=mtt.uniqueHit1; //Only one hit has highest score
			correctUniqueHit1+=mtt.correctUniqueHit1; //unique highest hit on answer site 
			correctMultiHit1+=mtt.correctMultiHit1;  //non-unique highest hit on answer site 
			correctLowHit1+=mtt.correctLowHit1;  //hit on answer site, but not highest scorer 
			noHit1+=mtt.noHit1;
			
			totalNumCorrect1+=mtt.totalNumCorrect1; //Skimmer only
			totalNumIncorrect1+=mtt.totalNumIncorrect1; //Skimmer only
			totalNumIncorrectPrior1+=mtt.totalNumIncorrectPrior1; //Skimmer only
			totalNumCapturedAllCorrect1+=mtt.totalNumCapturedAllCorrect1; //Skimmer only
			totalNumCapturedAllCorrectTop1+=mtt.totalNumCapturedAllCorrectTop1; //Skimmer only
			totalNumCapturedAllCorrectOnly1+=mtt.totalNumCapturedAllCorrectOnly1; //Skimmer only
			
			perfectMatch1+=mtt.perfectMatch1; //Highest slow score is max slow score
			semiperfectMatch1+=mtt.semiperfectMatch1; //A semiperfect mapping was found
			perfectMatchBases1+=mtt.perfectMatchBases1;
			semiperfectMatchBases1+=mtt.semiperfectMatchBases1;
			
			duplicateBestAlignment1+=mtt.ambiguousBestAlignment1;

			initialSiteSum1+=mtt.initialSiteSum1;
			postTrimSiteSum1+=mtt.postTrimSiteSum1;
			postRescueSiteSum1+=mtt.postRescueSiteSum1;
			siteSum1+=mtt.siteSum1;
			topSiteSum1+=mtt.topSiteSum1;
			
			AbstractIndex index=mtt.index();
			callsToScore+=index.callsToScore;
			callsToExtend+=index.callsToExtendScore;
			initialKeys+=index.initialKeys;
			initialKeyIterations+=index.initialKeyIterations;
			usedKeys+=index.usedKeys;
			usedKeyIterations+=index.usedKeyIterations;
			
			for(int j=0; j<index.hist_hits.length; j++){
				int x=Tools.min(hist_hits.length-1, j);
				hist_hits[x]+=index.hist_hits[j];
				hist_hits_score[x]+=index.hist_hits_score[j];
				hist_hits_extend[x]+=index.hist_hits_extend[j];
			}
			
			matchCountS1+=mtt.matchCountS1;
			matchCountI1+=mtt.matchCountI1;
			matchCountD1+=mtt.matchCountD1;
			matchCountM1+=mtt.matchCountM1;
			matchCountN1+=mtt.matchCountN1;

			mapped2+=mtt.mapped2;
			mappedRetained2+=mtt.mappedRetained2;
			rescuedP2+=mtt.rescuedP2;
			rescuedM2+=mtt.rescuedM2;
			lowQualityReadsDiscarded2+=mtt.lowQualityReadsDiscarded2;
			truePositiveP2+=mtt.truePositiveP2;
			truePositiveM2+=mtt.truePositiveM2;
			falsePositive2+=mtt.falsePositive2;
//			System.err.println("Adding "+mtt.falsePositive+" false positives -> "+falsePositive);
			totalCorrectSites2+=mtt.totalCorrectSites2;

			firstSiteCorrectP2+=mtt.firstSiteCorrectP2;
			firstSiteCorrectM2+=mtt.firstSiteCorrectM2;
			firstSiteIncorrect2+=mtt.firstSiteIncorrect2;
			firstSiteCorrectLoose2+=mtt.firstSiteCorrectLoose2;
			firstSiteIncorrectLoose2+=mtt.firstSiteIncorrectLoose2;
			firstSiteCorrectPaired2+=mtt.firstSiteCorrectPaired2;
			firstSiteCorrectSolo2+=mtt.firstSiteCorrectSolo2;
			firstSiteCorrectRescued2+=mtt.firstSiteCorrectRescued2;
			
			perfectHit2+=mtt.perfectHit2; //Highest score is max score
			perfectHitCount2+=mtt.perfectHitCount2;
			semiPerfectHitCount2+=mtt.semiPerfectHitCount2;
			uniqueHit2+=mtt.uniqueHit2; //Only one hit has highest score
			correctUniqueHit2+=mtt.correctUniqueHit2; //unique highest hit on answer site 
			correctMultiHit2+=mtt.correctMultiHit2;  //non-unique highest hit on answer site 
			correctLowHit2+=mtt.correctLowHit2;  //hit on answer site, but not highest scorer 
			noHit2+=mtt.noHit2;
			
			totalNumCorrect2+=mtt.totalNumCorrect2; //Skimmer only
			totalNumIncorrect2+=mtt.totalNumIncorrect2; //Skimmer only
			totalNumIncorrectPrior2+=mtt.totalNumIncorrectPrior2; //Skimmer only
			totalNumCapturedAllCorrect2+=mtt.totalNumCapturedAllCorrect2; //Skimmer only
			totalNumCapturedAllCorrectTop2+=mtt.totalNumCapturedAllCorrectTop2; //Skimmer only
			totalNumCapturedAllCorrectOnly2+=mtt.totalNumCapturedAllCorrectOnly2; //Skimmer only
			
			perfectMatch2+=mtt.perfectMatch2; //Highest slow score is max slow score
			semiperfectMatch2+=mtt.semiperfectMatch2; //A semiperfect mapping was found
			perfectMatchBases1+=mtt.perfectMatchBases1;
			semiperfectMatchBases1+=mtt.semiperfectMatchBases1;
			
			duplicateBestAlignment2+=mtt.ambiguousBestAlignment2;

			initialSiteSum2+=mtt.initialSiteSum2;
			postTrimSiteSum2+=mtt.postTrimSiteSum2;
			postRescueSiteSum2+=mtt.postRescueSiteSum2;
			siteSum2+=mtt.siteSum2;
			topSiteSum2+=mtt.topSiteSum2;
			
			matchCountS2+=mtt.matchCountS2;
			matchCountI2+=mtt.matchCountI2;
			matchCountD2+=mtt.matchCountD2;
			matchCountM2+=mtt.matchCountM2;
			matchCountN2+=mtt.matchCountN2;
			
		}
		maxReads=readsUsed1;
		if(syntheticReads>0){SYNTHETIC=true;}
		
		t.stop();
		long nanos=t.elapsed;
		
		if(verbose_stats>1){
			StringBuilder sb=new StringBuilder(1000);
			sb.append("\n\n###################\n#hits\tcount\tscore\textend\n");
			for(int i=0; i<hist_hits.length; i++){
				sb.append(i+"\t"+hist_hits[i]+"\t"+hist_hits_score[i]+"\t"+hist_hits_extend[i]+"\n");
			}
			try {
				ReadWrite.writeString(sb, "hist_hits.txt", true);
			} catch (Throwable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		final long basesUsed=(basesUsed1+basesUsed2);

		final double invTrials=1d/maxReads;
		final double invTrials100=100d/maxReads;
		double invSites100=100d/siteSum1;

		final double matedPercent=(numMated*invTrials100);
		final double badPairsPercent=(badPairs*invTrials100);
		final double innerLengthAvg=(innerLengthSum*1d/numMated);
		final double outerLengthAvg=(outerLengthSum*1d/numMated);
		final double insertSizeAvg=(insertSizeSum*1d/numMated);
		
		final double readsPerSecond=((readsUsed1+readsUsed2)*1000000000d)/nanos;
		final double fragsPerSecond=(keysUsed*1000000000d)/nanos;
		final double kiloBasesPerSecond=(basesUsed*1000000d)/nanos;
		
		double perfectHitPercent=(perfectHit1*invTrials100); //Highest score is max score
		double perfectMatchPercent=(perfectMatch1*invTrials100);
		double semiperfectMatchPercent=(semiperfectMatch1*invTrials100);
		
		double perfectHitCountPercent=perfectHitCount1*invSites100;
		double semiPerfectHitCountPercent=semiPerfectHitCount1*invSites100;
		
		double uniqueHitPercent=(uniqueHit1*invTrials100); //Only one hit has highest score
		double correctUniqueHitPercent=(correctUniqueHit1*invTrials100); //unique highest hit on answer site 
		double correctMultiHitPercent=(correctMultiHit1*invTrials100);  //non-unique highest hit on answer site 
		double correctLowHitPercent=(correctLowHit1*invTrials100);  //hit on answer site, but not highest scorer 
		double ambiguousFound=(duplicateBestAlignment1*invTrials100);
		double correctHighHitPercent=((correctMultiHit1+correctUniqueHit1)*invTrials100);
		double correctHitPercent=((correctLowHit1+correctMultiHit1+correctUniqueHit1)*invTrials100);

		double mappedB=(mapped1*invTrials100);
		double mappedRetainedB=(mappedRetained1*invTrials100);
		double rescuedPB=(rescuedP1*invTrials100);
		double rescuedMB=(rescuedM1*invTrials100);
		double falsePositiveB=(firstSiteIncorrect1*invTrials100);
		double falsePositiveLooseB=(firstSiteIncorrectLoose1*invTrials100);
		double truePositivePB=(firstSiteCorrectP1*invTrials100);
		double truePositiveMB=(firstSiteCorrectM1*invTrials100);
		double truePositiveStrict=((firstSiteCorrectP1+firstSiteCorrectM1)*invTrials100);
		double truePositiveLoose=(firstSiteCorrectLoose1*invTrials100);
		double snrStrict=10*Math.log10((firstSiteCorrectM1+firstSiteCorrectP1+0.1)/(firstSiteIncorrect1+0.1));
		double snrLoose=10*Math.log10((firstSiteCorrectLoose1+0.1)/(firstSiteIncorrectLoose1+0.1));
		double truePositivePMRatio=(truePositivePB/truePositiveMB);
		double truePositivePairedB=(firstSiteCorrectPaired1*100d/numMated);
		double truePositiveSoloB=(firstSiteCorrectSolo1*100d/(mappedRetained1-numMated));
		double truePositiveRescuedB=(firstSiteCorrectRescued1*100d/(rescuedP1+rescuedM1));
		double noHitPercent=(noHit1*invTrials100);
		
		long mappedReads, unambiguousReads;
		if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
			mappedReads=mappedRetained1+duplicateBestAlignment1;
			unambiguousReads=mappedRetained1;
		}else{
			mappedReads=mappedRetained1;
			unambiguousReads=mappedRetained1-duplicateBestAlignment1;
		}
		
		double avgNumCorrect=(SKIMMER ? totalNumCorrect1*invTrials : (totalCorrectSites1/(1d*(truePositiveP1+truePositiveM1))));
		double avgNumIncorrect=totalNumIncorrect1*invTrials; //Skimmer only
		double avgNumIncorrectPrior=totalNumIncorrectPrior1*invTrials; //Skimmer only

		double rateCapturedAllCorrect=totalNumCapturedAllCorrect1*invTrials100; //Skimmer only
		double rateCapturedAllTop=totalNumCapturedAllCorrectTop1*invTrials100; //Skimmer only
		double rateCapturedAllOnly=totalNumCapturedAllCorrectOnly1*invTrials100; //Skimmer only

		double avgCallsToScore=(callsToScore*invTrials);
		double avgCallsToExtendScore=(callsToExtend*invTrials);
		double avgInitialKeys=(initialKeys*1d/initialKeyIterations);
		double avgUsedKeys=(usedKeys*1d/usedKeyIterations);
		
		double avgInitialSites=(initialSiteSum1*invTrials);
		double avgPostTrimSites=(postTrimSiteSum1*invTrials);
		double avgPostRescueSites=(postRescueSiteSum1*invTrials);
		double avgSites=(siteSum1*invTrials);
		double avgPerfectSites=(perfectHitCount1*invTrials);
		double avgSemiPerfectSites=(semiPerfectHitCount1*invTrials);
		double avgTopSites=(topSiteSum1*invTrials);
		double lowQualityReadsDiscardedPercent=(lowQualityReadsDiscarded1*invTrials100);

		long matchErrors=matchCountS1+matchCountI1+matchCountD1;
		long baseLen=matchCountM1+matchCountI1+matchCountS1+matchCountN1;
		long matchLen=matchCountM1+matchCountI1+matchCountS1+matchCountN1+matchCountD1;
		long refLen=matchCountM1+matchCountS1+matchCountN1+matchCountD1;
		double errorRate=matchErrors*100d/matchLen;
		double matchRate=matchCountM1*100d/matchLen;//baseLen;
		double subRate=matchCountS1*100d/matchLen;//baseLen;
		double delRate=matchCountD1*100d/matchLen;
		double insRate=matchCountI1*100d/matchLen;//baseLen;
		double nRate=matchCountN1*100d/matchLen;//baseLen;
		
		if(SYNTHETIC && verbose_stats==-1){verbose_stats=Tools.max(verbose_stats,9);}

		sysout.println("Reads_Used"+DELIMITER+(readsUsed1+readsUsed2));
		sysout.println("Bases_Used"+DELIMITER+(basesUsed));
		sysout.println(String.format("Reads/sec"+DELIMITER+"%.2f", readsPerSecond));
		sysout.println(String.format("kBases/sec"+DELIMITER+"%.2f", kiloBasesPerSecond));
		double milf=msaIterationsLimited*invTrials;
		double milu=msaIterationsUnlimited*invTrials;
		if(verbose_stats>=1){sysout.println("MSA_iterations"+DELIMITER+String.format("%.2fL + %.2fU = %.2f", milf,milu,milf+milu));}
		
//		sysout.println();
//		sysout.println("\nRead 1 data:");
		
		sysout.println();
		
		if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
			double x=ambiguousFound+mappedRetainedB;
			sysout.println("R1_Mapped_Percent"+DELIMITER+padPercentMachine(x,4)+"%");
			sysout.println("R1_Unambiguous_Percent"+DELIMITER+padPercentMachine(mappedRetainedB,4)+"%");
			sysout.println("R1_Mapped_Reads"+DELIMITER+mappedReads);
			sysout.println("R1_Unambiguous_Reads"+DELIMITER+unambiguousReads);
		}else{			
			double x=mappedRetainedB-ambiguousFound;
			sysout.println("R1_Mapped_Percent"+DELIMITER+padPercentMachine(mappedRetainedB,4)+"%");
			sysout.println("R1_Unambiguous_Percent"+DELIMITER+padPercentMachine(x,4)+"%");
			sysout.println("R1_Mapped_Reads"+DELIMITER+mappedReads);
			sysout.println("R1_Unambiguous_Reads"+DELIMITER+unambiguousReads);
		}
		
		sysout.println();
		if(paired){
			sysout.println(String.format("Mated_Pairs"+DELIMITER+"%.4f%%", matedPercent));
			sysout.println(String.format("Bad_Pairs"+DELIMITER+"%.3f%%", badPairsPercent));
		}
		if(paired){
			sysout.println(String.format("R1_Rescued"+DELIMITER+"%.3f", rescuedPB+rescuedMB)+"%");
			sysout.println(String.format("Avg_Insert_Size"+DELIMITER+"%.2f", insertSizeAvg));
		}
		sysout.println();
		sysout.println(String.format("R1_Perfect_Best_Site"+DELIMITER+"%.4f", perfectMatchPercent)+"%");
		sysout.println(String.format("R1_Semiperfect_Site"+DELIMITER+"%.4f", semiperfectMatchPercent)+"%");
		sysout.println(String.format("R1_Ambiguous_Mapping"+DELIMITER+"%.4f", ambiguousFound)+"%");
//				+(REMOVE_DUPLICATE_BEST_ALIGNMENTS ? " (Removed)" : " (Kept)"));
		sysout.println(String.format("R1_Low_Quality_Discards"+DELIMITER+"%.4f", lowQualityReadsDiscardedPercent)+"%");
		
		if(MAKE_MATCH_STRING){
			sysout.println();
			sysout.println("R1_Match_Rate"+DELIMITER+padPercentMachine(matchRate,4)+"%");
			sysout.println("R1_Error_Rate"+DELIMITER+padPercentMachine(errorRate,4)+"%");
			sysout.println("R1_Sub_Rate"+DELIMITER+padPercentMachine(subRate,4)+"%");
			sysout.println("R1_Del_Rate"+DELIMITER+padPercentMachine(delRate,4)+"%");
			sysout.println("R1_Ins_Rate"+DELIMITER+padPercentMachine(insRate,4)+"%");
			sysout.println("R1_N_Rate"+DELIMITER+padPercentMachine(nRate,4)+"%");
			
			sysout.println("R1_Match_Count"+DELIMITER+matchCountM1);
			sysout.println("R1_Error_Count"+DELIMITER+matchErrors);
			sysout.println("R1_Sub_Count"+DELIMITER+matchCountS1);
			sysout.println("R1_Del_Count"+DELIMITER+matchCountD1);
			sysout.println("R1_Ins_Count"+DELIMITER+matchCountI1);
			sysout.println("R1_N_Count"+DELIMITER+matchCountN1);
		}
		
		if(paired){
			invSites100=100d/siteSum2;
			
			perfectHitPercent=perfectHit2*invTrials100; //Highest score is max score
			perfectMatchPercent=perfectMatch2*invTrials100;
			semiperfectMatchPercent=semiperfectMatch2*invTrials100;
			
			perfectHitCountPercent=perfectHitCount2*invSites100;
			semiPerfectHitCountPercent=semiPerfectHitCount2*invSites100;
			
			uniqueHitPercent=uniqueHit2*invTrials100; //Only one hit has highest score
			correctUniqueHitPercent=correctUniqueHit2*invTrials100; //unique highest hit on answer site 
			correctMultiHitPercent=correctMultiHit2*invTrials100;  //non-unique highest hit on answer site 
			correctLowHitPercent=correctLowHit2*invTrials100;  //hit on answer site, but not highest scorer 
			ambiguousFound=(duplicateBestAlignment2*invTrials100);
			correctHighHitPercent=(correctMultiHit2+correctUniqueHit2)*invTrials100;
			correctHitPercent=(correctLowHit2+correctMultiHit2+correctUniqueHit2)*invTrials100;

			mappedB=(mapped2*invTrials100);
			mappedRetainedB=(mappedRetained2*invTrials100);
			rescuedPB=(rescuedP2*invTrials100);
			rescuedMB=(rescuedM2*invTrials100);
			falsePositiveB=(firstSiteIncorrect2*invTrials100);
			falsePositiveLooseB=(firstSiteIncorrectLoose2*invTrials100);
			truePositivePB=(firstSiteCorrectP2*invTrials100);
			truePositiveMB=(firstSiteCorrectM2*invTrials100);
			truePositiveStrict=((firstSiteCorrectP2+firstSiteCorrectM2)*invTrials100);
			truePositiveLoose=(firstSiteCorrectLoose2*invTrials100);
			snrStrict=10*Math.log10((firstSiteCorrectM2+firstSiteCorrectP2+0.1)/(firstSiteIncorrect2+0.1));
			snrLoose=10*Math.log10((firstSiteCorrectLoose2+0.1)/(firstSiteIncorrectLoose2+0.1));
			truePositivePMRatio=(truePositivePB/truePositiveMB);
			truePositivePairedB=(firstSiteCorrectPaired2*100d/numMated);
			truePositiveSoloB=(firstSiteCorrectSolo2*100d/(mappedRetained2-numMated));
			truePositiveRescuedB=(firstSiteCorrectRescued2*100d/(rescuedP2+rescuedM2));
			avgNumCorrect=(totalCorrectSites2/(1d*(truePositiveP2+truePositiveM2)));
			noHitPercent=noHit2*invTrials100;
			
			avgNumCorrect=(SKIMMER ? totalNumCorrect2*invTrials : (totalCorrectSites2/(1d*(truePositiveP2+truePositiveM2))));
			avgNumIncorrect=totalNumIncorrect1*invTrials; //Skimmer only
			avgNumIncorrectPrior=totalNumIncorrectPrior1*invTrials; //Skimmer only

			rateCapturedAllCorrect=totalNumCapturedAllCorrect2*invTrials100; //Skimmer only
			rateCapturedAllTop=totalNumCapturedAllCorrectTop2*invTrials100; //Skimmer only
			rateCapturedAllOnly=totalNumCapturedAllCorrectOnly2*invTrials100; //Skimmer only
			
			if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
				mappedReads=mappedRetained2+duplicateBestAlignment2;
				unambiguousReads=mappedRetained2;
			}else{
				mappedReads=mappedRetained2;
				unambiguousReads=mappedRetained2-duplicateBestAlignment2;
			}

			avgInitialSites=initialSiteSum2*invTrials;
			avgPostTrimSites=postTrimSiteSum2*invTrials;
			avgPostRescueSites=postRescueSiteSum2*invTrials;
			avgSites=siteSum2*invTrials;
			avgPerfectSites=(perfectHitCount1*invTrials);
			avgSemiPerfectSites=(semiPerfectHitCount1*invTrials);
			avgTopSites=topSiteSum2*invTrials;
			lowQualityReadsDiscardedPercent=lowQualityReadsDiscarded2*invTrials100;

			matchErrors=matchCountS2+matchCountI2+matchCountD2;
			baseLen=matchCountM2+matchCountI2+matchCountS2+matchCountN2;
			matchLen=matchCountM2+matchCountI2+matchCountS2+matchCountN2+matchCountD2;
			refLen=matchCountM2+matchCountS2+matchCountN2+matchCountD2;
			errorRate=matchErrors*100d/matchLen;
			matchRate=matchCountM2*100d/matchLen;//baseLen;
			subRate=matchCountS2*100d/matchLen;//baseLen;
			delRate=matchCountD2*100d/matchLen;
			insRate=matchCountI2*100d/matchLen;//baseLen;
			nRate=matchCountN2*100d/matchLen;//baseLen;
			
//			sysout.println("\n\nRead 2 data:");
			sysout.println();
//			sysout.println(String.format("perfectHit"+DELIMITER+"%.2f", perfectHitPercent)+"%");
//			sysout.println(String.format("uniqueHit"+DELIMITER+"%.2f", uniqueHitPercent)+"%");
//			sysout.println(String.format("correctUniqueHit"+DELIMITER+"%.2f", correctUniqueHitPercent)+"%");
//			sysout.println(String.format("correctMultiHit"+DELIMITER+"%.2f", correctMultiHitPercent)+"%");
//			sysout.println(String.format("correctHighHit"+DELIMITER+"%.2f", correctHighHitPercent)+"%");
//			sysout.println(String.format("correctHit"+DELIMITER+"%.2f", correctHitPercent)+"%");
			
			//sysout.println(String.format("mapped"+DELIMITER+(mappedB<10?" ":"")+"%.3f", mappedB)+"%");
			if(REMOVE_DUPLICATE_BEST_ALIGNMENTS){
				double x=ambiguousFound+mappedRetainedB;
				sysout.println("R2_Mapped_Percent"+DELIMITER+padPercentMachine(x,4)+"%");
				sysout.println("R2_Unambiguous_Percent"+DELIMITER+padPercentMachine(mappedRetainedB,4)+"%");
				sysout.println("R2_Mapped_Reads"+DELIMITER+mappedReads);
				sysout.println("R2_Unambiguous_Reads"+DELIMITER+unambiguousReads);
			}else{			
				double x=mappedRetainedB-ambiguousFound;
				sysout.println("R2_Mapped_Percent"+DELIMITER+padPercentMachine(mappedRetainedB,4)+"%");
				sysout.println("R2_Unambiguous_Percent"+DELIMITER+padPercentMachine(x,4)+"%");
				sysout.println("R2_Mapped_Reads"+DELIMITER+mappedReads);
				sysout.println("R2_Unambiguous_Reads"+DELIMITER+unambiguousReads);
			}
			sysout.println();
			if(paired){
				sysout.println(String.format("R2_Rescued"+DELIMITER+"%.3f", rescuedPB+rescuedMB)+"%");
			}
			sysout.println();
			sysout.println(String.format("R2_Perfect_Best_Site"+DELIMITER+"%.4f", perfectMatchPercent)+"%");
			sysout.println(String.format("R2_Semiperfect_Site"+DELIMITER+"%.4f", semiperfectMatchPercent)+"%");
			sysout.println(String.format("R2_Ambiguous_Mapping"+DELIMITER+"%.4f", ambiguousFound)+"%");
								//(REMOVE_DUPLICATE_BEST_ALIGNMENTS ? "(Removed)" : "(Kept)"));
			sysout.println(String.format("R2_Low_Quality_Discards"+DELIMITER+"%.4f", lowQualityReadsDiscardedPercent)+"%");
			
			if(MAKE_MATCH_STRING){
				sysout.println();
				sysout.println("R2_Match_Rate"+DELIMITER+padPercentMachine(matchRate,4)+"%");
				sysout.println("R2_Error_Rate"+DELIMITER+padPercentMachine(errorRate,4)+"%");
				sysout.println("R2_Sub_Rate"+DELIMITER+padPercentMachine(subRate,4)+"%");
				sysout.println("R2_Del_Rate"+DELIMITER+padPercentMachine(delRate,4)+"%");
				sysout.println("R2_Ins_Rate"+DELIMITER+padPercentMachine(insRate,4)+"%");
				sysout.println("R2_N_Rate"+DELIMITER+padPercentMachine(nRate,4)+"%");
				
				sysout.println("R2_Match_Count"+DELIMITER+matchCountM2);
				sysout.println("R2_Error_Count"+DELIMITER+matchErrors);
				sysout.println("R2_Sub_Count"+DELIMITER+matchCountS2);
				sysout.println("R2_Del_Count"+DELIMITER+matchCountD2);
				sysout.println("R2_Ins_Count"+DELIMITER+matchCountI2);
				sysout.println("R2_N_Count"+DELIMITER+matchCountN2);
			}
		}
		
		if(BBSplitter.TRACK_SCAF_STATS){
			BBSplitter.printCounts(BBSplitter.SCAF_STATS_FILE, BBSplitter.scafCountTable, true, readsUsed1+readsUsed2);
		}
		
		if(BBSplitter.TRACK_SET_STATS){
			BBSplitter.printCounts(BBSplitter.SET_STATS_FILE, BBSplitter.setCountTable, true, readsUsed1+readsUsed2);
		}
		
		errorState|=ReadStats.writeAll(paired);
		
		assert(!CALC_STATISTICS || truePositiveP1+truePositiveM1+falsePositive1+noHit1+lowQualityReadsDiscarded1==maxReads) : 
			"\nThe number of reads out does not add up to the number of reads in.\nThis may indicate that a mapping thread crashed.\n"+
			truePositiveP1+"+"+truePositiveM1+"+"+falsePositive1+"+"+noHit1+"+"+lowQualityReadsDiscarded1+" = "+
			(truePositiveP1+truePositiveM1+falsePositive1+noHit1+lowQualityReadsDiscarded1)+" != "+maxReads;
		if(!SKIMMER){
			assert(!CALC_STATISTICS || truePositiveP1+truePositiveM1==correctLowHit1+correctMultiHit1+correctUniqueHit1);
		}else{
			assert(!CALC_STATISTICS || truePositiveP1+truePositiveM1==correctLowHit1+correctUniqueHit1);
		}
	}
	
	static final void printSettings0(int k, int maxindel, float minratio){
		if(MACHINE_OUTPUT){
			sysout.println("Genome"+DELIMITER+Data.GENOME_BUILD);
			sysout.println("Key_Length"+DELIMITER+k);
			sysout.println("Max_Indel"+DELIMITER+maxindel);
			sysout.println("Minimum_Score_Ratio"+DELIMITER+minratio);
			sysout.println("Mapping_Mode"+DELIMITER+(PERFECTMODE ? "perfect" : SEMIPERFECTMODE ? "semiperfect" : "normal"));
		}else{
			sysout.println("Genome:                \t"+Data.GENOME_BUILD);
			sysout.println("Key Length:            \t"+k);
			sysout.println("Max Indel:             \t"+maxindel);
			sysout.println("Minimum Score Ratio:  \t"+minratio);
			sysout.println("Mapping Mode:         \t"+(PERFECTMODE ? "perfect" : SEMIPERFECTMODE ? "semiperfect" : "normal"));
		}
	}
	
	
	static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}
	
	protected static void clearStatics(){
		maxReads=-1;
//		readsUsed=0;
//		readsUsed2=0;
//		lowQualityReadsDiscarded1=0;
//		lowQualityReadsDiscarded2=0;
//		lowQualityBasesDiscarded1=0;
//		lowQualityBasesDiscarded2=0;
		
		outputBaseName="readsOut_"+(System.nanoTime()&0x1FFFF);
		outFile=null;
		outFile2=null;
		outFileM=null;
		outFileM2=null;
		outFileU=null;
		outFileU2=null;
		outFileB=null;
		outFileB2=null;
		ArrayList<String> blacklist=null;
		
		errorState=false;
	}
	
	/* ------------ Non-static fields ----------- */
	

	ConcurrentReadInputStream cris;
	ConcurrentReadOutputStream rosA=null, rosM=null, rosU=null, rosB=null;
	
	float fractionGenomeToExclude=-1;
	int maxIndel1=-1;
	int maxIndel2=-1;
	int minApproxHits=-1;
	int expectedSites=-1;
	int ambigMode=AMBIG_BEST;
//	int ambigMode2=AMBIG_BEST;
	boolean fast=false;
	boolean slow=false;
	boolean vslow=false;
	boolean verbose=false;
	boolean rcompMate=false;
	boolean outputSitesOnly=false;
	long targetGenomeSize=-1;
	int ziplevel=-1;
	int build=1;
	String reference=null;
	int keylen=13;
	int idmodulo=1;
	float samplerate=1f;
	double minid=-1;
	long sampleseed=1;
	boolean ambiguousRandom=false, ambiguousAll=false;
	boolean forceanalyze=false;
	private boolean gunzip=false;
	private boolean gzip=false;
	private boolean pigz=false;
	private boolean unpigz=false;
	boolean setxs=false, setintron=false;
	String bamscript=null;
	String in1=null, in2=null, qfin1=null, qfin2=null;
	String qfout=null, qfout2=null, qfoutM=null, qfoutM2=null, qfoutU=null, qfoutU2=null, qfoutB=null, qfoutB2=null;
	
	/** Scores below the (max possible alignment score)*(MINIMUM_ALIGNMENT_SCORE_RATIO) will be discarded.
	 * Default: 0.4 ~ 0.5 for clean data against raw PacBio data.
	 * Very sensitive!  A value of 0.2 will potentially produce many false positives. */
	float MINIMUM_ALIGNMENT_SCORE_RATIO;

	float keyDensity;//Normal key density
	float maxKeyDensity; //For situations where some of the read is too low quality, this is the max for the rest of the read. 
	float minKeyDensity;
	int maxDesiredKeys; //Don't go above this number of keys except to maintain minKeyDensity.
	
	/** Additional ref bases on each end of site mapping location in alignment window.
	 * If there are no insertions or deletions, 0 is fine. */
	int SLOW_ALIGN_PADDING;
	int SLOW_RESCUE_PADDING;
	int TIP_SEARCH_DIST;
	
	/** Class name of MSA to use */
	String MSA_TYPE;
	int MAX_SITESCORES_TO_PRINT;
	boolean PRINT_SECONDARY_ALIGNMENTS;
	
	/* ------------ Coverage ----------- */
	
	CoveragePileup pileup;
	String coverageStats=null, coverageBinned=null, coverageBase=null, coverageHist=null;
	int coverageMinScaf=0;
	boolean cov32bit=false;
	boolean covBitset=false;
	boolean covSetbs=false;
	boolean covArrays=true;
	boolean covNzo=false;
	boolean covTwocolumn=false;
//	boolean covSecondary=true;
	boolean covKsb=true;
	boolean covStranded=false;
	boolean covStartOnly=false;
	int covBinSize=1000;
	
	
	/* ------------ Static fields ----------- */

	static final int AMBIG_BEST=0;
	static final int AMBIG_TOSS=1;
	static final int AMBIG_RANDOM=2;
	static final int AMBIG_ALL=3;
	
	static int THRESH=0; //Threshold for calculating true positives on synthetic data, or something. 
	
	static int synthReadlen=150;

	static int maxInsLen=30; //Default 40
	static int maxSubLen=30; //Default 40
	static int maxDelLen=40; //Default 8000
	
	static byte minQuality=3;
	static byte midQuality=23;
	static byte maxQuality=35;
	
	static int maxSnps=4;//4;
	static int maxInss=3;//2;
	static int maxDels=3;
	static int maxSubs=3;//2;
	
	static float baseSnpRate=0.50f;
	static float baseInsRate=0.30f;
	static float baseDelRate=0.30f;
	static float baseSubRate=0.30f;//0.3f;
	static float PERFECT_READ_RATIO=0.0f;//0.2f;//0.8f
	
	//Extra work for rare cases in human only.
	static boolean SAVE_AMBIGUOUS_XY=false;
	
	static boolean colorspace=false;
	
	static boolean translateToBaseSpace=false; //Translate (colorspace) reads before outputting them
	

	static boolean TRIM_LIST=true; //Increases speed many times; reduces accuracy a bit

	static boolean PAIRED_RANDOM_READS=false;
	static boolean REQUIRE_CORRECT_STRANDS_PAIRS=true;
	static boolean SAME_STRAND_PAIRS=false;
	static boolean KILL_BAD_PAIRS=false;
	
	static boolean INDEX_LOADED=false;
	static final boolean SLOW_ALIGN=true; //Do a more accurate scoring pass with MSA
	static boolean MAKE_MATCH_STRING=SLOW_ALIGN;
	
	/** Rescue paired reads by searching near mate */
	static boolean RESCUE=true;
	
	/** Generally should be set to false unless SLOW_ALIGN==true */
	static boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;

	/** Forbid alignments with indels longer than MAX_INDEL */
	static boolean STRICT_MAX_INDEL=false;
	/** Don't allow reads to map to their origin location in the reference. Useful for self-correcting reads. */
	static boolean FORBID_SELF_MAPPING=false;
	/** Only allow perfect and semiperfect mappings */
	static boolean SEMIPERFECTMODE=false;
	/** Only allow perfect mappings */
	static boolean PERFECTMODE=false;
	/** Only allow sites with at least this many contiguous matches */
	static int KFILTER=-1;
	/** Only allow sites with identity of at least this */
	static float IDFILTER=0f;
	
	/** Quality-trim left side of read before mapping */
	static boolean qtrimLeft=false;
	/** Quality-trim right side of read before mapping */
	static boolean qtrimRight=false;
	/** Restore read to untrimmed state after mapping (and destroy match string) */
	static boolean untrim=false;
	/** Trim bases with quality less than or equal to this value */
	static byte TRIM_QUALITY=6;
	/** Don't trim reads to be shorter than this */
	static int minTrimLength=60;
	/** Produce local alignments instead of global alignments */
	static boolean LOCAL_ALIGN=false;
	
	public static int minChrom=1;
	public static int maxChrom=Integer.MAX_VALUE;

	static long maxReads=-1;
	
	protected static boolean CALC_STATISTICS=true;

	static boolean QUICK_MATCH_STRINGS=false;
	static boolean OUTPUT_READS=false;
	static boolean DONT_OUTPUT_UNMAPPED_READS=false;
	static boolean DONT_OUTPUT_BLACKLISTED_READS=false;
	
	static boolean OUTPUT_ORDERED_READS=false;
	static boolean DOUBLE_PRINT_ERROR_RATE=false;
	static boolean PRINT_UNMAPPED_COUNT=false;
	
	static String outputBaseName="readsOut_"+(System.nanoTime()&0x1FFFF);
	static String outFile=null;//outputBaseName+"_1.txt";
	static String outFile2=null;//outputBaseName+"_2.txt";
	static String outFileM=null;//outputBaseName+"_mapped_1.txt";
	static String outFileM2=null;//outputBaseName+"_mapped_2.txt";
	static String outFileU=null;//outputBaseName+"_unmapped_1.txt";
	static String outFileU2=null;//outputBaseName+"_unmapped_2.txt";
	static String outFileB=null;//outputBaseName+"_blacklist_1.txt";
	static String outFileB2=null;//outputBaseName+"_blacklist_2.txt";
	static ArrayList<String> blacklist=null;
	static ArrayList<String> splitterOutputs=null;

	static boolean useRandomReads=false;
	static int sequentialOverlap=5;
	static boolean sequentialStrandAlt=false;

	static boolean overwrite=false;
	static boolean append=false;
	static boolean SYNTHETIC=false;
	static boolean ERROR_ON_NO_OUTPUT=false;
	static boolean MACHINE_OUTPUT=false;
	static boolean USE_MODULO=false;
	final static String DELIMITER="=";
	
	static PrintStream sysout=System.err;
	static boolean SYSIN=false;
	static int verbose_stats=0;
	static boolean waitForMemoryClear=false;
	
	public static boolean errorState=false;
	
}
