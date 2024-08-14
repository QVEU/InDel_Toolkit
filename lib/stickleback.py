#!/usr/bin/env python
'''
stickleback.py
v0.2
  ><```º>
Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit
3/31/23
USAGE: python stickleback.py pathto/input.sam query templateFasta [] []   ( [arg] = optional argument)
Designed to map inserted sequences in a template sequence (plasmid or transcript).
v0.2 version notes:
    Long inserts created an issue with mapping, and highlight a issue specific to nanopore, that a good match for the insert may be flanked by bad sequence or adapters.
    Fixed by limiting the flanking sequence to 25 for all hits (regardless of the size of the insert)
'''

##### Imports #####
import pandas as pd
import Levenshtein
import time
import numpy as np
import sys
from multiprocessing import Pool

##### Functions #####
def Initialize(args):
    '''
    Function: Initialize()
        Reads in and filters data and initializes a bunch of parameters for the mapping.

    Arguments:
        "seq_query" is a *TUPLE* of the template sequence and the query. Need to be sent as a package for multi-threading with Pool.
    Note: replaces position if alternative minimum sites are found. Producing a 3' bias? Need to consider improvements here.
    '''
    #Input SAM
    print("\n\n----------------=============-----------------")
    print("--==--==--==--==   ><```º>   ==--==--==--==--=")
    print("==--==--==--==-- stickleback --==--==--==--==-")
    print("----------------=============-----------------\n\n")


    #Check for Args
    if len(args)<3:
        print("USAGE: python stickleback.py pathto/input.sam query templateFasta [minLength] [maxLength]   ( [arg] = optional argument)")
        exit()

    #Query
    query=(args[2].upper())
    print("Query Length: {}".format(len(query)))

    #Template File
    template=args[3]
    with open(template,'r') as fasta:
        templateSeq="".join([i.strip() for i in fasta][1:])
    print("Template Length: {}".format(len(templateSeq)))
    
    Dtemp=blockDist((templateSeq.upper(),query.upper()))[0] #compute the distance to the template (without insert) to determine sensitivity.
    print("Template-Query Distance: {}".format(Dtemp))

    #Choose Read Size
    if len(args)>4:
        minL=int(args[4]) # read length minimum
        maxL=int(args[5]) # read length maximum
    else:
        minL=len(query)
        maxL=len(templateSeq)

    cutoff=(int(Dtemp)-2)
    print("Mapping Reads of Size {} to {} with a cutoff of {}".format(minL, maxL, cutoff))

    #input SAM
    inputData=sys.argv[1]
    print("\nLoading SAM: {}".format(inputData))
    pdSam=loadSAM(inputData,minL,maxL)
    outfile=inputData.replace(".sam",'_stickleback.csv')
    print("\nOutfile: {}".format(outfile))

    return(str(query),str(templateSeq), pdSam, int(Dtemp), int(cutoff), outfile)

def loadSAM(inputData, minL, maxL):
    '''
    Function: loadSAM()
        Reads in SAM file.
        Since this is the slow step, I broke it out to work on later.

    Arguments:
        "inputData" is a string of the path to the input SAM file.
    Note: quite slow, how can we speed up? pd.read seems to have issues with the headers on the file.
    '''
    with open(inputData,"r") as IF:
        samInput=[lin.split("\t")[0:11] for lin in IF if lin.startswith("@") is False]#Filter out headers with @

    names=["read","FLAG","template","pos","mapq","cigar","Rnext","Pnext","Tlen","seq","Qscore"]#SAM columns

    #types=[str,int,str,int,int,str,str,str,int,str,int]#SAM columns
    #typeDict=dict(zip(names,types))

    size=len(samInput)
    print("Total Candidate Reads: {}".format(size))
    pdSam=pd.DataFrame(samInput,columns=names)
    #print(pdSam)
    #pdSam.astype(typeDict)
    pdSam['length']=pdSam.seq.apply(lambda s: len(s))
    pdSam=pdSam[(pdSam.length>minL)&(pdSam.length<maxL)&(pdSam.template!="*")]
    return(pdSam)

def blockDist(seq_query):
    contextWindow=25
    '''
    ~ new in v0.1 ~
    Function: blockDist()
        unpacks the sequence query and template, makes chunks of query length from template and then computes the distance for each kmer.
        Uses numpy vectorize to speed up the computation (formerly for-loop)

    Arguments:
        `seq_query` is a tuple of read and query string, passed this way for 'pool'
    '''
    #Unpack Values
    seq=seq_query[0]
    query=seq_query[1]
    #print(seq)
    #Prepare Query
    quL=len(query)
    minD=quL #set minimum distance to length of query

    # Build vectorized
    BlockDistance = lambda Block: Levenshtein.distance(Block,query)
    vectorBlockDist = np.vectorize(BlockDistance,otypes=[int])

    # Make blocks of read sequence
    blocks=np.array([seq[i:(i+quL)] for i in range(len(seq)-quL)])
    bDist = vectorBlockDist(blocks)
    #compute minimum and position of minimum
    minD=min(bDist.astype(int))
    minPos=np.argmin(bDist.astype(int))

    matchlist=minPos
    matchseq=blocks[minPos]
    contextSeq="|".join([seq[max(0,minPos-contextWindow):minPos],seq[min(len(seq),minPos+quL):min(len(seq),(minPos+(quL+contextWindow)))]])#made a change here 3/30/23 PTD

    return(minD,minPos,matchseq,contextSeq,matchlist)

def poolBlocks(query,pdSam,nthreads=14): #uses multithreading to compute the matches (Edited to 14 threads 4/27/23 WB)
    with Pool(nthreads) as p:
        print("\n1. Computing minimum distance hit position for {} reads.".format(len(pdSam)))
        quL=len(query)
        print("Query: {}".format(query))
        pdSam['minD'],   pdSam['minPos'],   pdSam['matchseq'],    pdSam['context'],   pdSam['matchlist']   = zip(*p.map(blockDist, [(c.upper(),query.upper()) for c in pdSam.seq]))

        print("2. Mapping minimum distance site on template sequence...")
        pdSam['minD_v'], pdSam['minPos_v'], pdSam['matchseq_v'] , pdSam['context_v'], pdSam['matchlist_v'] = zip(*p.map(blockDist, [(templateSeq.upper(),c.replace("\\|","").upper()) for c in pdSam.context]))
        print("done.")
        pdSam['insPos_v']=pdSam["minPos_v"]+26 #EDITED 4/22/23 PTD. SAME AS contextWindow EDITEDWB:+26to fix python indexing
    return(pdSam)

if __name__=="__main__":
    t0=(time.time())
    query, templateSeq, pdSam, Dtemp, cutoff, outfileName = Initialize(sys.argv)
    pdSam = poolBlocks(query,pdSam)
    print("\nMapped hits in {} reads.".format(np.sum([int(i)<int(cutoff) for i in pdSam.minD])))
    pdSam[(pdSam.minPos>0)&(pdSam.minD<cutoff)].to_csv(outfileName)
    t1=(time.time())
    print("Done in {} minutes.".format((t1-t0)/60))
    print("Wrote {}.".format(outfileName))
