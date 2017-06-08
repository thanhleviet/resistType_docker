import sys
from Bio import SeqIO

records = SeqIO.index(sys.argv[1], "fasta")
outputFile=sys.argv[1].split(".fa")[0] + "_padded.fa"
with open(outputFile, "w") as f:
    for record in sorted(records):
        sequence = str(records[record].seq)
        n=100
        newSequence = 'N'*100 + sequence + 'N'*100
        f.write(">{0}\n{1}\n".format(record, "\n".join([newSequence[i:i+n] for i in range(0, len(newSequence), n)])))
                

