'''
Hang Phan
7/7/2015
This script read a nucleotide fasta file containing resistance gene database (many information, including resistance profile) 
and create several files including forward DNA, protein sequence, renaming genes so that it is computer friendly. 
'''
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
promoters = ['P3', 'P4', 'P5', 'Pa/Pb', 'Pc/Pd', 'ampC_promoter', 'ampC promoter']
foAllGene = open("resistDB/ResistanceGeneSeqs.fasta", "w")   
foAA = open("resistDB/ResistanceGeneSeqs_aa.fasta", "w")
foFW = open("resistDB/ResistanceGeneSeqs_fw.fasta", "w")
foOtherDNA = open("resistDB/ResistanceGeneSeqs_otherDNA.fasta", "w")
mapGeneName = open("resistDB/ResistanceGeneSeqs_mapGeneName.txt", "w")
foResistProfile=open("resistDB/ResistanceProfiles.txt", "w")
geneList=[]
f=open(sys.argv[1], "r")
for line in f:
    cols=line.strip().split("\t")
    seqid = cols[0]
    dnaSeq = cols[11].strip().rstrip().replace("-", "").replace(" ", "")
    protSeq = cols[14].strip().rstrip()
    if seqid == "gene_name":
        foResistProfile.write(line)
        continue

    geneName = seqid.replace("(", "").replace(")", "").replace("\'", "_prime").replace("/", "_").replace(",", "_").replace(" ", "").replace(" ", "_")
    newLine = "\t".join([geneName] + cols[1:]) + "\n"
    foResistProfile.write(newLine)

    if seqid == "SHV-100":
        foAllGene.write(">{0}\n{1}\n".format(seqid, dnaSeq))
        mapGeneName.write("{0}\t{1}\n".format(seqid, geneName))    
        continue
    if geneName not in geneList:
        geneList.append(geneName)
    else:
        print "duplicated gene", geneName
        print "add this one as a separate entry"
        geneName += "_1"

    mapGeneName.write("{0}\t{1}\n".format(seqid, geneName))    
    
    if dnaSeq.strip().rstrip() == "" and protSeq != "":
        foAA.write(">{0}\n{1}\n".format(geneName, protSeq))
    elif dnaSeq.strip().rstrip() != "":
        foAllGene.write(">{0}\n{1}\n".format(geneName, dnaSeq))
        seqLen = len(dnaSeq)
        codonStart = seqLen % 3
        if seqid in promoters:
            foFW.write(">{0}\n{1}\n".format(geneName, dnaSeq))
            foOtherDNA.write(">{0}\n{1}\n".format(geneName, dnaSeq))
            continue
        expectedProtLen = seqLen/3
        ok = False
        for codonStart in range(3):
            protSeq = Seq(str(dnaSeq)[codonStart:], generic_dna).translate(table = 11, to_stop=True)
            if len(protSeq) > expectedProtLen - 4:
                foFW.write(">{0}\n{1}\n".format(geneName, dnaSeq))
                foAA.write(">{0}\n{1}\n".format(geneName, protSeq))
                ok=True
                break
        if not ok:
            for codonStart in range(3):
                protSeq = Seq(str(Seq(dnaSeq).reverse_complement())[codonStart:], generic_dna).translate(table =11, to_stop = True)
                if len(protSeq) > expectedProtLen - 4:
                    foFW.write(">{0}\n{1}\n".format(geneName, Seq(dnaSeq).reverse_complement()))
                    foAA.write(">{0}\n{1}\n".format(geneName, protSeq)) 
                    ok = True
                    break
            if not ok:
                print seqid, "protein sequence not detected properly. Ignore the generation of protein sequence."
                foFW.write(">{0}\n{1}\n".format(geneName, dnaSeq))
                foOtherDNA.write(">{0}\n{1}\n".format(geneName, dnaSeq))

foAA.close()
foFW.close()
foOtherDNA.close()
f.close()
foAllGene.close()
foResistProfile.close()
