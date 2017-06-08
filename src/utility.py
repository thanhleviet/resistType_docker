import os
from Bio import SeqIO
import numpy as np
import math


def readConfig(configFile):
    '''
    get information of the configurations (where softwares can be found)
    '''
    

                    
def checkAssemblySize(contigFile):
    '''
    Given an assembly, calculate the total size of the assembly, see if it's contaminated (or mixed sample), single sample or metagenomic samples
    '''
    if not os.path.exists(contigFile):
        return [0,0]

    records = SeqIO.index(contigFile, "fasta")
    lens = np.array([len(records[record].seq) for record in records])
    sumLen= np.sum(lens)
    minSize=2.5*1000000
    maxSize=7.0*1000000
    nSamples=0

    if sumLen < minSize:
        nSamples=0
    if sumLen >= minSize and sumLen <= maxSize:
        nSamples=1
    elif sumLen >= maxSize * 1.6 and sumLen <=maxSize*2.2:
        nSamples = 2
    elif sumLen > maxSize * 2.2:
        nSamples = sumLen * 2 / ( maxSize + minSize)
    
    return [nSamples, sumLen]

def getAssemblyMeanCoverage(contigFile):
    '''
    Given an assembly, get the mean coverage of the assembly from the contig name description
    '''
    records= SeqIO.index(contigFile, "fasta")
    covs = []
    lens = []
    allLens=0.0
    for record in records:
        descriptions =[record]
        lenSeq = len(records[record].seq)
        allLens += lenSeq
        cov = getCovContig(record, records)
        if cov and cov > 1 and lenSeq > 250 and cov < 1000: 
            covs.append(cov)
            lens.append(len(records[record].seq))

    if len(covs) ==0:
        return [0,0,0]

    meanCov, stdCov = weighted_avg_and_std(np.array(covs), np.array(lens))
    newCovs=[]
    newLens=[]
    for i in range(0, len(covs)):
        if covs[i] > meanCov - 2*stdCov and covs[i] < meanCov + 2* stdCov:
            newCovs.append(covs[i])
            newLens.append(lens[i])
    newMeanCov, newStdCov= weighted_avg_and_std(np.array(newCovs), np.array(newLens))
    print "Mean Coverage ", newMeanCov, newStdCov 
    return [newMeanCov, newStdCov, np.sum(np.array(newLens))/allLens]
        
def getCovContig(record, records):
    descriptions = [record] + records[record].description.split()
    cov = None
    for i in descriptions:
        if "_cov_" in i:
            cols = i.split("_")
            cov =float( cols[cols.index("cov")+1])
            break
    return cov

def estimateContigCopyNumber(record, records,  meanCov, stdCov):
    thisCov = getCovContig(record, records)

    nonRoundedRatio=thisCov/meanCov
    return round(nonRoundedRatio,2)
    ratio =round(nonRoundedRatio, 0)

    
    if thisCov > ratio * meanCov - stdCov *nonRoundedRatio and thisCov < ratio * meanCov + stdCov *nonRoundedRatio:
        return ratio
    elif nonRoundedRatio < ratio:
        if ratio ==1:
            return 1
        else:
            return ratio -1
    else:
        return ratio
    
        
def getSTDContigs(records, meanCov, stdCov):
    '''
    get large contigs >20kb that have coverage within the mean +- std values
    '''
    outContigs = []
    totalLen = 0

    for record in records:
        
        thisCov = getCovContig(record, records)
        lenRecord = len(records[record].seq)
        if thisCov > meanCov - 2*stdCov and thisCov < meanCov + 2*stdCov  and lenRecord > 10000:
            outContigs.append([record, lenRecord])
        
    return outContigs
    
    
     
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))




def mergeIntervals(intervals):
    '''
    merge a list of intervals into none overlapping ones
    '''

    if not intervals:
        return []
    intervals = sorted(intervals, key = lambda x: x[0])
    result = []
    (a, b) = intervals[0]
    for (x, y) in intervals[1:]:
        if x <= b:
            b = max(b, y)
        else:
            result.append((a, b))
            (a, b) = (x, y)
    result.append((a, b))
    return result
