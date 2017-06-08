import sys, os
from Bio import SeqIO

clusters = {}
clusterFile = "resistDB/temp90.fa.clstr"
fastaFile = "resistDB/ResistanceGeneSeqs_fw.fasta"
records = SeqIO.index(fastaFile, "fasta")
chromGenes=['gyrA_ecol', 'gyrB_ecol', 'gyrA_kpne', 'gyrB_kpne', 'parC_ecol', 'parC_kpne', 'parE_ecol', 'parE_kpne']
genes, repGene = [], None
with open(clusterFile, "r") as f:
    for line in f:
        if line.startswith(">"):
            if len(genes) > 0:
                clusters[repGene] = genes
                genes = []
        else:
            geneName = line.strip().split()[2].rstrip(".")[1:]

            genes.append(geneName)
            if line.strip().endswith("*"):
                repGene = geneName

    clusters[repGene] = genes
for repGene in clusters.keys():
    if len(clusters[repGene]) >=1:
        f = open("resistDB/clusterFas/{0}.fa".format(repGene), 'w')
        for gene in clusters[repGene]:
            f.write(">{0}\n{1}\n".format(gene, records[gene].seq))
        f.close()
