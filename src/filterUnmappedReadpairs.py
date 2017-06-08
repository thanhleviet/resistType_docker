import sys, os, gzip
'''
input from stdin is in sam formated sorted by name. This script tries to removed unmapped read pairs only and put the remaining read pairs into 2 fastq files.
'''
fi = sys.stdin
prevReadName = ""

firstInPair  = []
secondInPair = []

fq1 = gzip.open(sys.argv[1], "w")
fq2 = gzip.open(sys.argv[2], "w")



for line in fi:
    if line.startswith("@"):
        #print line.strip()
        pass
    else:
        cols = line.strip().split()
        readName = cols[0]
        flag = int(cols[1])
        
        if flag & 4 and flag & 8:
            continue
        
        if prevReadName == "" or readName  != prevReadName:
            if firstInPair != [] and secondInPair!=[]:
                fq1.write("@{0}\n{1}\n+{0}\n{2}\n".format(firstInPair[0], firstInPair[1], firstInPair[2]))
                fq2.write("@{0}\n{1}\n+{0}\n{2}\n".format(secondInPair[0], secondInPair[1], secondInPair[2]))
                pass            
            firstInPair = []
            secondInPair = []
            prevReadName = readName

        if flag&64 :#read is first in pair
            firstInPair=[cols[0], cols[9], cols[10]] #get name, sequence and quality line
        if flag & 128: #read is second in pair
            secondInPair =[cols[0], cols[9], cols[10]] #get name, sequence and quality line

if firstInPair != [] and secondInPair!=[]:
    fq1.write("@{0}\n{1}\n+{0}\n{2}\n".format(firstInPair[0], firstInPair[1], firstInPair[2]))
    fq2.write("@{0}\n{1}\n+{0}\n{2}\n".format(secondInPair[0], secondInPair[1], secondInPair[2]))
    
fq1.close()
fq2.close()

'''
for line in fi:
    if line.startswith("@"):
        print line.strip()
    else:
        cols = line.strip().split()
        readName = cols[0]
        flag = int(cols[1])

        if flag & 4 and flag & 8:
            continue

        if prevReadName == "" or readName  != prevReadName:
            if stackedLines!= []:
                for item in stackedLines:
                    print item.strip()
                pass
            stackedLines = [line]            
            prevReadName = readName
        else:
            stackedLines.append(line)
for item in stackedLines:
    print item.strip()
'''

