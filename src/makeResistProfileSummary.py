'''
Author
Hang Phan - hang.phan@ndm.ox.ac.uk
11 Jan  2016
Given a list of samples, look in to the folders containing resistPred.txt and get the resistance gene profile, this version contain the copy number report 
'''
from __future__ import division

import sys, os, pysam, gzip, logging, subprocess, uuid, shutil
import operator
import logging.handlers
import time, datetime
import os.path
from optparse import  OptionParser

_baseDir = '/well/bag/hangphan/resistType/'

# Set up logging
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger('Log')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)
MODE1, MODE2, MODE3, MODE4 = 0,1,2,3
class makeResistGeneProfiles(object):
    def __init__(self, args):
        self.virulenceFactorList=['afaE1_SaT040','afaE3','astA_ecol_42','bmaE_IH11165','cdtB_IHE3034','chuA_CFT073','chuS_CFT073','chuT_CFT073','chuU_CFT073','chuW_CFT073','chuX_CFT073','chuY_CFT073','clbA_IHE3034','clbB_IHE3034','clbN_IHE3034','clbQ_IHE3034','clpG','cnf1','cvaC_pAPECO2ColV','draA_IH11128','draB_IH11128','draC_IH11128','draD_IH11128','draE_IH11128','draP_IH11128','entA_CFT073','entB_CFT073','entC_CFT073','entD_CFT073','entE_CFT073','entF_CFT073','F17-like_fimbrial_Chaperone_ecol_536','F17-like_fimbrial_subuniT_ecol_536','F17-like_fimbrial_usher_ecol_536','F17-like_fimbrial_adhesin_subunit_ecol_536','feoB_CFT073','fepA_CFT073','fepB_CFT073','fepC_CFT073','fepD_CFT073','fepE_CFT073','fepG_CFT073','fimA_CFT073','fimB_CFT073','fimC_CFT073','fimD_CFT073','fimE_CFT073','fimF_CFT073','fimG_CFT073','fimH_CFT073','fimI','focA_CFT073','focC_CFT073','focD_CFT073','focF_CFT073','focG_CFT073','focH_CFT073','sfaD_CFT073','fyuA_APEC01','gafD','H7_fliC_CFT073','hlyA_CFT073','hlyB_CFT073','hlyC_CFT073','hlyD_CFT073','hlyF_pETN48','hra_AY040859.1_Hek','ibeA_VF0237','ibeB_RS218','ibeC_VF0237','iha_O157_H7','ireA_APEC01','iroN_CFT073','irp-2_APEC01','iss_AF042279.1','iucA_CFT073','iucB_CFT073','iucC_CFT073','iucD_CFT073','iutA_CFT073','kfiC','kpsD_CFT073','kpsM_CFT073','kpsM_K15_ecol_536','kpsT','malX_CFT073','nfaE_O83_K1_H4','ompA_NMEC58C','ompT_CFT073','papA_CFT073','papB_CFT073','papC_CFT073','papD_CT073','papE_CFT073','papF_CFT073','papG_CFT073','papH_CFT073','papI_CFT073','papJ_CFT073','papK_CFT073','pic_CFT073','rfc_K12_MG1655','sat_ecolCFT073','sfaA_VF0235','sfaB_ecol_536','sfaC_ecol_536','sfaD_ecol_536','sfaE_ecol_536','sfaF_ecol_536','sfaG_ecol_536','sfaH_ecol_536','sfaS_ecol_536','tcpC_CFT073','traJ_VF0241','traT_pAPECO2ColV','tsh_CFT073','uspA_CFT073','vat_APEC01','yfcV_CFT073']
        self.sampleFile=args[0]
        self.modePred=args[1]
        self.modeGene=args[2]
        self.header=None
        self.drugDict={}
        self.drugList=[]
        self.samples=[]
        self.outputFN1 = self.sampleFile.split(".txt")[0] + "_phenotypePred.csv"
        self.outputFN2 = self.sampleFile.split(".txt")[0] + "_genotypePred.csv"
        self.geneList=['CTX-M', 'TEM', 'TEM-promoter', 'OXA', 'SHV', 'aac', 'aad', 'aph', 'ant', 'KPC', 'CMY', 'NDM', 'VIM', 'tet', 'gyrA']
        self.temPromoter=['P3', 'Pa_Pb', 'Pc_Pd', 'P4', 'P5', 'Pa/Pb', 'Pc/Pd']
    def makeHeaderPred(self): #make the header line for the phenotypic prediction
        header = ""
        
        if self.modePred == MODE1:
            with open("{0}/resistDB/drugnamesInResistDB.txt".format(_baseDir), "r") as f:
                for idx, line in enumerate(f):
                    self.drugDict[line.strip()] = 1 + idx
                    self.drugList.append( line.strip())
            self.header = ",".join(["sampleID"] + self.drugList)

        if self.modePred == MODE3:
            with open("{0}/resistDB/jacResults_header.txt".format(_baseDir), "r") as f:
                for idx, line in enumerate(f):
                    self.drugDict[line.strip()] = 1 + idx
                    self.drugList.append( line.strip())
            self.header = ",".join(["sampleID"] + self.drugList)
        if self.modePred == MODE4:
            with open("{0}/resistDB/drugList2.txt".format(_baseDir), "r") as f:
                for idx, line in enumerate(f):
                    self.drugDict[line.strip()] = 1 + idx
                    self.drugList.append( line.strip())
            self.header = ",".join(["sampleID"] + self.drugList)
    
        if self.modePred == MODE2:
            with open("{0}/resistDB/headerForHICF.csv".format(_baseDir), "r") as f:
                header = f.readlines()[0]
                cols = header.strip().split(",")
                for idx, item in enumerate(cols):
                    if ( "Routine lab" in item or "Phoenix" in item ) and "SIR" in item:
                        drugName = item.replace("Routine lab ", "").replace("Phoenix ", "").split()[0]
                        if drugName in drugList + ["amoxicillin", "ampicillin", "ertapenem", "meropenem"]:
                            self.drugDict[drugName] = idx
            self.header = header
        with open(self.outputFN1, 'w') as f:
            f.write("{0}\n".format( self.header))

        return None

    def makeHeaderGene(self):
        if self.modeGene==0:
            header= ['sampleID'] + self.geneList + ['Others']
            
        else:
            header = ['sampleID', 'matchState', 'pident', 'geneName', 'resistanceProfile']

        with open(self.outputFN2, 'w') as f:
            f.write("{0}\n".format(",".join(header)))

    def run(self):
    
        with open(self.sampleFile, "r") as f:
            for line in f:
                self.samples.append(line.strip().split()[0])
        self.makeHeaderPred()
        self.makeHeaderGene()
        for sample in self.samples:
            exactList, inexactList, resistanceProfile = self.readResult(sample)
            if exactList ==None and inexactList==None and resistanceProfile ==None:
                continue
            self.writeGenePred(sample, exactList, inexactList)
            self.writePhenoPred(sample, resistanceProfile)
        
    def readResult(self, sample):
        resultFile = "resistType/{0}/resistancePred.txt".format(sample)
        if not os.path.exists(resultFile):
            logger.info("{0}'s resistance prediction not found \n".format(sample))
            return None, None, None
        exactList={}
        inexactList={}
        print resultFile
        with open(resultFile, "r") as f:

            lines = f.readlines()
            n=0
            exactMatch = 0
            inexactMatch = 0
            while n < len(lines):
                line = lines[n]
                n+=1
                if line.startswith("#"):
                    continue
                if "List of " in line:
                    continue
                if "Exact gene match" in line:
                    exactMatch = 1
                    continue
                if exactMatch:
                    if "Inexact" in line:
                        exactMatch = 0
                        inexactMatch = 1
                        continue
                    else:
                        cols = line.strip().rstrip().split(",")
                        geneName = cols[0].split(":")[0]
                        if geneName not in exactList:
                            exactList[geneName]= [line.strip()]
                        else:
                            exactList[geneName].append(line.strip())
                        
                if inexactMatch:
                    if "* Resistance prediction" in line:
                        break
                    cols = line.strip().rstrip().split(",")
                    if float(cols[1].split(":")[0].split("|")[0]) < 80:
                        continue
                                                        
                    inexactList[cols[0].split(":")[0]]=[line.strip()]
                if "* Resistance prediction" in line:
                    break

        resistanceProfile = {}
        resistanceTerms=["R", "( R )", "R(MIC>8)", "low level R"]
        sureResistanceTerms=["R", "R(MIC>8)"]
        for line in lines[n:]:
            cols = line.strip().replace(", possible", "-possible").split(",")
            #print line
            if len(cols) <=2:
                continue
            #print cols
            if ":" not in cols[2]:
                continue
            for x in cols[2:]:
                
                x = x.rstrip().strip()
                y = x.split(":")

                if y[0] not in resistanceProfile:
                    resistanceProfile[y[0]] = y[1]
                else:
                    if resistanceProfile[y[0] ] in  ['S', 'no literature']:
                        resistanceProfile[y[0]] = y[1]
                    if (resistanceProfile[y[0]] == '( R )' or "?" in resistanceProfile[y[0]]) and (y[1] in sureResistanceTerms):
                        resistanceProfile[y[0]] == y[1]

        return exactList, inexactList, resistanceProfile

            
    def writeGenePred(self, sample, exactList, inexactList):
        foG=  open(self.outputFN2, "a")
        nCols=100        
        if self.modeGene ==0:  #report a selected set of genes in separate columns, then aggregate the remaining into one single column, separated by |
            outLine = ['']*(len(self.geneList) + 2)
            outCols=[]*(len(self.geneList) + 2)
            for idx in range(len(self.geneList)+ 2):
                outCols.append([])
            outLine[0] =  sample
            outCols[0] = [sample]

            for item in exactList:
                inList=0
                if item in self.virulenceFactorList:
                    continue
                for idx, gene in enumerate(self.geneList):
                    if item.startswith(gene):
                        outCols[idx + 1].append(item)
                        inList=1
                        break
                if item in self.temPromoter:
                    idx = self.geneList.index("TEM-promoter")
                    outCols[idx +1].append(item )
                    inList=1
                if inList==0:
                    outCols[-1].append(item)
            for item in inexactList:
                if item in self.virulenceFactorList:
                    continue

                inList =0
                pident = inexactList[item][0].split(",")[1].split(":")[0].split("|")[0]
                #if float(pident) == 100:
                #    pident = "98.0"
                for idx, gene in enumerate(self.geneList):
                    if item.startswith(gene) and float(pident) >80:
                        outCols[idx + 1].append(item + ":" + pident)
                        inList=1
                        break
                if item in self.temPromoter:
                    idx = self.geneList.index("TEM-promoter")
                    outCols[idx +1].append(item + ":" + pident)
                    inList=1
                if inList==0:
                    outCols[-1].append(item + ":"+ pident)

            for idx in range(1, len(outCols)):
                outLine[idx] = "|".join(outCols[idx])
            foG.write("{0}\n".format(",".join(outLine)))

        else: #write one gene per line

            for item in sorted(exactList.keys()):
                if item in self.virulenceFactorList:
                    continue

                lines= exactList[item]
                for line in lines:
                    cols = line.strip().rstrip().split(",")
                    outCols = [""]*nCols
                    outCols[0] = sample
                    outCols[1] = "exactMatch"
                    outCols[2] = "100.0"
                    outCols[3] = cols[0].split(":")[0]
                    if len(cols[0].split(":")) ==1:
                        print sample, line
                    outCols[4] = cols[0].split(":")[1]

                    for idx in range(1, len(cols)):
                        outCols[idx +4] = cols[idx].rstrip().strip()
                    
                    foG.write("{0}\n".format(",".join(outCols)))
            for item in sorted(inexactList.keys()):
                if item in self.virulenceFactorList:
                    continue

                lines= inexactList[item]
                for line in lines:
                    cols = line.strip().rstrip().split(",")
                    outCols = [""]*nCols
                    outCols[0] = sample
                    outCols[1] = "inexactMatch"
                    outCols[2] = cols[1].strip()
                    outCols[3]= cols[0].split(":")[0]
                    outCols[4] = cols[0].split(":")[1]
                    for idx in range(2, len(cols)):
                        outCols[idx +3] = cols[idx].rstrip().strip()
                    foG.write("{0}\n".format(",".join(outCols)))
            foG.close()
        return
    def writePhenoPred(self, sample, resistanceProfile):
        '''
        write the resistance profile prediction to output file, using several options
        '''

        drugDict=self.drugDict
        nCols=len(drugDict.keys())+1
        line = ["S"] * nCols
        line[0] = sample
        if self.modePred==MODE1: # all antibiotics mode
            for x in resistanceProfile:
                if x == "virulence factor":
                    continue
                if x == "ampicillin":
                    line[self.drugDict["amoxicillin_ampicillin"]] = resistanceProfile[x]
                else:
                    line[self.drugDict[x]] = resistanceProfile[x]

        if self.modePred == MODE2: # only drugs in HICF dataset
            for x in resistanceProfile:
                if x in drugDict:
                    line[drugDict[x]] = resistanceProfile[x]
                if x == "amoxicillin_ampicillin":
                    line[drugDict["amoxicillin"]] = resistanceProfile[x]
                    line[drugDict["ampicillin"]] = resistanceProfile[x]
                if x == "carbapenems":
                    line[drugDict["ertapenem"]] = resistanceProfile[x]
                    line[drugDict["meropenem"]] = resistanceProfile[x]
                if x == "aminoglycosides and quinolones":
                    line[drugDict["aminoglycosides"]] = resistanceProfile[x]
                    line[drugDict["quinolones"]] = resistanceProfile[x]

        if self.modePred == MODE3: # only drugs in JAC dataset
            for x in resistanceProfile:
                if x in drugDict:
                    line[drugDict[x]] = resistanceProfile[x]
                if x == "amoxicillin_ampicillin":
                    line[drugDict["ampicillin"]] = resistanceProfile[x]
                if x == "carbapenems":
                    line[drugDict["ertapenem"]] = resistanceProfile[x]
                    line[drugDict["meropenem"]] = resistanceProfile[x]

        if self.modePred == MODE4: # only drugs in BloodCulture dataset
            for x in resistanceProfile:
                if x in drugDict:
                    line[drugDict[x]] = resistanceProfile[x]

        for idx, item in enumerate(line):
            if "does not confer resistance" in item:
                line[idx] = "S" 
                continue
            if "consider" in item:
                line[idx] = "S"
                continue
            if "no literature" in item or " may cause" in item:
                line[idx] = "NS" #not sure

        with open(self.outputFN1, 'a') as fo:
            fo.write("{0}\n".format(",".join(line)))
        return


if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-s", "--sampleFile", dest="sampleFile", type="string", default=None, help="File containing sample names")
    parser.add_option("-m", "--predMode", dest="outputPredictionMode", type="int", default=0, help="Mode of output for the resistance prediction \n\t\t\t0: all drug mode \n\t\t\t1:HICF comparison mode \n\t\t\t2: JAC set comparison mode \n\t\t\t3: BloodCultureSample mode")
    parser.add_option("-g", "--geneMode", dest="outputGeneMode", type="int", default=0, help="Mode of output for the resistance gene\n\t\t\t0: report selected genes mode, one row per sample\n\t\t\t1: report all genes mode, a sample cover several rows, with associated resistance mechanisms")
    (opts, args) = parser.parse_args()
    
    makeResistGeneProfileModule= makeResistGeneProfiles([opts.sampleFile, opts.outputPredictionMode, opts.outputGeneMode])
    makeResistGeneProfileModule.run()
