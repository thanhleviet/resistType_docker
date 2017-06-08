#!/usr/bin/env python
'''
MLST typing for a couple of species

Author
Hang Phan - hang.phan@ndm.ox.ac.uk
23 Apr 2016
Last updated 07/06/2017

If the assembly is available, run from the contigs, otherwise look for fastq files an run assembly of only MLST gene locus from scratch using spades.

'''
from __future__ import division
import sys, os, pysam, gzip, logging, subprocess, uuid, shutil
import operator
import logging.handlers
import time, datetime
import os.path
from subprocess import Popen, PIPE
from optparse import  OptionParser
import numpy as np
from pysam import Fastafile
from Bio import SeqIO


#Directory where src, resistDB, bin are sitting
_baseDir="/".join(os.path.dirname(os.path.realpath(sys.argv[0])).split("/")[:-1])

# Set up logging
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger('Log')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)

class mlstType(object):
    def __init__(self, args):
        self.args = args
        self.sampleid = args[0]
        self.mlstDir= _baseDir + "/mlstdb/"
        
        self.reducedMLSTAlleles = self.mlstDir + "temp90_padded.fa"
        self.alleleDB = self.mlstDir + "allAlleles.fa"
        self.stProfileFile = self.mlstDir + "allSTs.tsv"
        
        self.contigFile = args[1]
        self.outDir = args[2]
        self.fq1 = args[3]
        self.fq2 = args[4]

        if self.outDir == None:
            self.outDir = "mlst/"
        self.prefix = self.outDir + self.sampleid
        cmdLine = "mkdir -p " + self.outDir
        os.system(cmdLine)

        self.miniContigFile = "{0}.mlst.contigs.fa".format(self.prefix)
        if not self.contigFile:
            self.contigFile=self.miniContigFile
            
        self.blastFile = "{0}.outputBlast.txt".format(self.prefix)
        self.mlstOutFile = "{0}.mlstResults.txt".format(self.prefix)
        self.matchedAlleles = []
        self.closestMatchedAlleles=[]
        self.STProfiles={}
        self.AlleleSTMap={}
        
        self.bwaPath = _baseDir + "/bin/bwa"
        self.spadesPath= "spades.py"
        self.blastprog = _baseDir + "/bin/blastn"

    def runMLST(self):
        self.readSTProfiles()
        if not os.path.exists(self.contigFile):
            self.runMapAndSpades()
            self.contigFile = self.miniContigFile
        
        self.runBlast()
        self.getMatchedAlleles()
        self.makeMLSTPrediction()
        return
    def makeMLSTPrediction(self):
        '''
        make prediction of mlst from matched alleles 
        '''
        groupAllele = [i.split("-")[0] for i in self.matchedAlleles]
        import collections
        counter = collections.Counter(groupAllele)
        mainSpecies = []
        speciesPrefix = ""
        f = open(self.mlstOutFile, "w")
        
        for (group, count) in counter.most_common(2):
            
            alleleTypeA = [i for i in self.matchedAlleles if i.startswith(group)]
            alleleGroup = ["-".join(i.split("-")[0:2]) for i in alleleTypeA] 
            alleleCounter = collections.Counter(alleleGroup)

            if len(alleleCounter) > 4:
                mainSpecies.append(group)
                
            lenMostCommonAlleleGroup = alleleCounter.most_common(1)[0][1]

            if len(alleleCounter) ==7 and len(alleleTypeA) ==7:

                hashIdx = ".".join(sorted(alleleTypeA))
                if hashIdx in self.AlleleSTMap:
                    f.write("{0},{1},{2}\n".format(self.sampleid, self.AlleleSTMap[hashIdx], hashIdx.replace(".", ",")))
                else:
                    f.write("{0},{1}-newST,{2}\n".format(self.sampleid,  hashIdx.split("-")[0], hashIdx.replace(".", ",")))

            else:
                thisSet = set(alleleTypeA)
                hasST=0
                bestMatch=[0, []]
                if len(alleleCounter) <4:
                    continue
                for st in self.STProfiles:
                    stSet = self.STProfiles[st]
                    if stSet.issubset(thisSet):
                        f.write("{0},{1},{2}\n".format(self.sampleid, st, ",".join(sorted(list(stSet)))))
                        hasST =1
                        break
                    overlap = stSet.intersection(thisSet)
                    if len(overlap) > bestMatch[0]:
                        bestMatch[0] = len(overlap)
                        bestMatch[1].append(st)
                            
                if hasST==0:
                    if lenMostCommonAlleleGroup >1:
                        f.write("{0},{1}-mixedST,{2}\n".format(self.sampleid, alleleTypeA[0].split("-")[0], ",".join(sorted(alleleTypeA))))
                    else:
                        f.write("{0},{1}-newST_newAllele,{2}\n".format(self.sampleid, alleleTypeA[0].split("-")[0],",".join(sorted(alleleTypeA))))

                    #f.write("#bestMatches (with {0} loci matching): (check matchedAlleles.txt file for sequence of close matched  alleles)\n".format(bestMatch[0]))
                    #for item in bestMatch[1]:
                    #    f.write("#{0},{1},{2}\n".format(self.sampleid, item, ",".join(sorted(self.STProfiles[item]))))
                        
        f.close()
        return
        
    def runMapAndSpades(self):
        '''
        map reads to reduced allele file and use those reads to assemble into contigs
        '''
        fqs = [self.fq1]

        if self.fq2 != None:
            fqs.append(self.fq2)
        filterScript= _baseDir + "/src/filterUnmappedReadpairs.py"
        self.miniFq1 = "{0}.mlst_reads1.fq.gz".format(self.prefix)
        self.miniFq2 = "{0}.mlst_reads2.fq.gz".format(self.prefix)
        #cmdLine = "{5} mem {0} {1} | python {2} {3} {4}".format( self.reducedMLSTAlleles, " ".join(fqs),filterScript, self.miniFq1, self.miniFq2, self.bwaPath)
        cmdLine = "bowtie2 --fast -x {0}  -1 {1} -2 {2}|python {3} {4} {5} ".format(self.reducedMLSTAlleles, self.fq1, self.fq2 , filterScript, self.miniFq1, self.miniFq2)
        print cmdLine
        os.system(cmdLine)
        cmdLine = "{3} -1 {0} -2 {1} -o tmpDir/{2} --careful".format(self.miniFq1, self.miniFq2, self.sampleid, self.spadesPath)
        os.system(cmdLine)
        os.system("cp tmpDir/{0}/contigs.fasta {1}".format(self.sampleid, self.miniContigFile))
        os.system("rm -rf tmpDir/{0}".format(self.sampleid))
        return

    def readSTProfiles(self):
        '''
        Read the sequence type profile of all species from a big file into a dictionary
        '''
        for line in open(self.stProfileFile):
            cols = line.strip().split()
            st = cols[0]
            alleles = sorted(cols[1:8])
            alleleString= ".".join(alleles)
            self.AlleleSTMap[alleleString] = st
            self.STProfiles[st] = set(alleles)
        return

    def runBlast(self):
        '''
        run blast of contig file against the database of all alleles 
        '''
        cmdLine = "{3} -query {0} -db {1} -word_size 17  -evalue 0.001  -gapopen 5 -gapextend 2 -culling_limit 1 -outfmt '6 qseqid sseqid pident length slen qstart qend sstart send qseq sseq'|awk '$3>95 && $4/$5>0.95' > {2}".format(self.contigFile, self.alleleDB, self.blastFile, self.blastprog)
        logger.info(cmdLine)
        os.system(cmdLine)
        return

    def getMatchedAlleles(self):
        '''
        read blast result from blastFile and get matched allele or it's closest hit
        '''
        self.matchedAlleles = []
        self.closestMatchedAlleles = []
        contig_allele = {}
        for line in open(self.blastFile):
            cols = line.strip().split()
            qseqid, sseqid = cols[:2]
            pident, length, slen, qstart, qend, sstart, send = map(float, cols[2:9])
            qseq, sseq = cols[9:11]

            if send < sstart:
                temp = send
                send = sstart
                sstart = temp
            lenRatio = length / (send - sstart +1)
            if lenRatio <0:  lenRatio = -lenRatio
            if lenRatio >1: lenRatio = 1/lenRatio
            pident = pident  * lenRatio
            newVector = [sseqid, pident, qstart, qend, qseq, sseq]
            if qseqid not in contig_allele:
                contig_allele[qseqid] = [newVector]
            else:
                #check overlap
                isOverlap = 0
                for idx, item in enumerate(contig_allele[qseqid]):
                    overlapSize = min(qend, item[3]) - max(qstart, item[2])
                    if overlapSize > slen *0.8 :
                        isOverlap = 1
                        if newVector[1] > item[1]:
                            contig_allele[qseqid][idx] = newVector
                if isOverlap ==0:
                    contig_allele[qseqid].append(newVector)
        self.matchedAlleleFile =  "{0}.matchedAlleles.txt".format(self.prefix)
        f = open(self.matchedAlleleFile, 'w')
        for item in contig_allele.keys():
            matches = contig_allele[item]
            for match in matches:
                f.write("{0},{1},{2},{3}\n".format(self.sampleid, match[0], match[1], match[-2], match[-1]))
                if match[1] ==100.0:
                    self.matchedAlleles.append(match[0])
                else:
                    self.closestMatchedAlleles.append(match)
        f.close()
        return
                    



if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-s", "--sampleName", dest="sampleName", type="string", default=None, help="Sample name")
    parser.add_option("-c", "--contigFile", dest="contigFile", type="string", default=None, help="link to bam file")
    parser.add_option("-1", "--fastqFile1", dest="fq1", type="string", default=None, help="link to fastqFile1")
    parser.add_option("-2", "--fastqFile2", dest="fq2", type="string", default=None, help="link to fastqFile2")
    parser.add_option("-o", "--outDir", dest = "outDir", type="string", default = None, help="directory of output file")
    (opts, args) = parser.parse_args()
    mlstModule= mlstType([opts.sampleName, opts.contigFile, opts.outDir, opts.fq1, opts.fq2])
    mlstModule.runMLST()

