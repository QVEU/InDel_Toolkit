#!/usr/bin/env python
# coding: utf-8
'''
stickleback.py
v0.2
  ><```º>
Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit
3/31/23
#usage: python DelMapper_v0.1.py <inputfile (str)> <delSize (int)> <outputfile (str)>

Designed to map deletions of defined sizes on a template sequence (plasmid or transcript).

v0.2 version notes:
    Long inserts created an issue with mapping, and highlight a issue specific to nanopore, that a good match for the insert may be flanked by bad sequence or adapters.
    Fixed by limiting the flanking sequence to 25 for all hits (regardless of the size of the insert)

GitHub push (8/14/24):
- cleaned up comments and checked out some repairs suggested by reviewer.
- Generated readmes and example data to accompany github repo. 

'''
 

TEMPLATE="ttaaaacagcctgtgggttgcacccacccacagggcccactgggcgctagcactctggtactgaggtacctttgtgcgcctgtttttactccccttcccccgaagtaacttagaagctgtaaatcaacgatcaatagcaggtgtggcacaccagtcataccttgatcaagcacttctgtttccccggactgagtatcaataggctgctcgcgcggctgaaggagaaaacgttcgttacccgaccaactacttcgagaagcttagtaccaccatgaacgaggcagggtgtttcgctcagcacaaccccagtgtagatcaggctgatgagtcactgcaacccccatgggcgaccatggcagtggctgcgttggcggcctgcccatggagaaatccatgggacgctctaattctgacatggtgtgaagagcctattgagctagctggtagtcctccggcccctgaatgcggctaatcctaactgcggagcacatgctcacaaaccagtgggtggtgtgtcgtaacgggcaactctgcagcggaaccgactactttgggtgtccgtgtttccttttattcctatattggctgcttatggtgacaatcaaagagttgttaccatatagctattggattggccatccggtgtgcaacagggcaattgtttacctatttattggttttgtaccattatcactgaagtctgtgatcactctcaaattcattttgaccctcaacacaatcaaacatgggctcacaagtgtccacacaacgctccggttcacacgaaaactctaactcagctaccgagggttccactataaactatactaccattaattactataaagattcctatgccgccacagcaggtaagcagagccttaagcaggacccagacaagtttgcaaatcctgtcaaagacatcttcactgaaatggcagcgccattaaaatctccatctgctgaggcatgtggttacagcgatcgggtggcacaattaactattggcaattctaccatcactacgcaagaagcagcaaacatcatagttggctatggtgagtggccttcctactgttcggactctgatgctactgcagtggacaaaccaacgcgcccagatgtttcggtgaataggttttacacattggacacaaaattgtgggagaaatcatccaaggggtggtactggaaattcccggatgtgttaactgaaaccggggtctttggtcaaaatgcacagttccactacctctatcggtcagggttctgcattcacgtgcagtgcaatgctagtaagttccaccaaggagcactcctagtcgctgtcctcccagagtatgtcattgggacagtggcaggtggcacagggacggaggatagccaccccccttataagcagactcaacccggtgctgatggcttcgaattgcaacacccgtacgtgcttgatgctggcattccaatatcacaattaacagtgtgcccacatcagtggattaatttgaggaccaacaattgtgccacaataatagtgccgtacataaacgcactaccctttgattctgccttgaaccattgtaactttggtctgctggttgtgcctattagcccgttagattatgaccaaggtgcgacgccagtgatccccattactatcactttggccccaatgtgttctgaatttgcaggccttagacaagcagttacgcaagggtttcctactgagctgaaacctggcacaaaccaatttttaaccactgacgatggcgtctcagcacccattctgccaaactttcaccccaccccgtgtatccatatacccggtgaagttagaaacttgctagagctatgccaggtggagaccattttagaggtcaacaatgtacctacgaatgccactagcttaatggagagactgcgcttcccggtctcagctcaagccgggaaaggtgagctatgtgcagtgttcagagctgaccctggacgaagtgggccatggcagtccaccttgttgggccagttgtgcgggtactacacccaatggtcaggatcactggaagtcaccttcatgttcaccgggtcctttatggctaccggcaagatgctcatagcatacacaccaccaggaggccccttacccaaggaccgggcgaccgccatgttgggcacgcacgtcatctgggactttgggctgcaatcgtctgtcacccttgtaataccatggatcagcaacactcattacagagcgcacgctcgagatggtgtgttcgactactacactacaggtttggttagcatatggtaccagacgaattatgtggttccaattggggcacccaatacagcctatataatagcattggcggcagcccagaagaacttcaccatgaagttgtgtaaggatgctagtgatatcctacagacaggcactatccagggagatagggtggcagatgtgattgagagttctataggggacagtgtgagcagagccctcacccgagctctaccggcacctaccggccaagacacacaggtaagcagccaccgattagatactggtaaagttccagcactccaagccgctgaaattggagcatcatcaaatgctagtgatgagagtatgattgagacacggtgtgttcttaattcacatagtacagctgagaccactcttgatagcttcttcagcagagcaggattagttggagagatagacctccctcttgaaggcacaaccaacccgaatgggtacgcaaactgggacatagacataacaggttacgcgcaaatgcgtagaaaggtggagctgttcacctacatgcgttttgacgcagagttcacctttgttgcatgcacccctaccgggcaagttgtcccgcaattgctccaatacatgtttgtaccacccggagcccccaagccagactccagagaatctctcgcatggcaaactgccactaatccctcagtttttgtgaagctgtcagaccccccagcacaggtttctgttccattcatgtcacctgcgagcgcctatcaatggttttatgacgggtatcccacattcggtgaacacaaacaggagaaagaccttgaatacggggcatgcccaaacaacatgatgggtacgttctcagtgcggactgtaggcacctcgaagtccaagtacccattggtgatcaggatttacatgaggatgaagcacgtcagggcgtggatacctcgcccaatgcgtaaccagaactatctattcaaagccaacccaaattatgctggtaattttattaaaccaactggtgccagtcgcacagcaatcaccaccctcgggaaatttggacagcagtccggagctatctacgtgggcaactttagagtggttaaccgccatcttgctactcataatgactgggcaaaccttgtttgggaagacagctcccgcgacttgctcgtatcatctaccactgctcaaggttgtgacacgattgctcgttgcaattgccagacaggagtgtattattgtaactcaatgagaaaacactatccggtcagtttctcgaaacccagtttgatcttcgtggaggccagcgagtattatccagctagataccagtcacatctcatgcttgcagtgggtcattcggaaccaggggattgcggtggcattcttagatgccaacatggcgtcgtagggatagtttccaccgggggaaacggcctggtggggttcgccgatgtgagggatcttctgtggttggatgatgaagccatggagcagggcgtgtctgattacattaaagggcttggagatgcttttggcatggggtttacagacgcagtgtcaagagaagttgaagcactgaaaagtcacttgatcggctcagagggtgccgtggagaagattctaaagaacttagttaaactcatctctgcgctcgtcatcgtcatcaggagtgattatgacatggtcacattgacggcaacacttgccctgatcgggtgccacgggagcccttgggcctgggttaagtcgaagacagcatcaatcttgggcataccgatggctcagaagcagagtgcctcttggttaaagaagttcaacgatgcggcgagtgccgcgaaggggcttgagtggatctccaacaaaatcagtaaatttatcgattggctcaaggagaaaatcataccggctgctaaagagaaagtcgagtttctaaacaatctaaagcaactccccttattggagaaccaaatttctaatctcgaacagtcagcagcttcgcaggaggaccttgaggcgatgtttggcaacgtgtcttatctggcccacttctgccgcaaattccaacccctctatgccacggaagcaaagagggtgtacgccctagaaaagagaatgaataattacatgcagttcaagagcaaacaccgtattgaacctgtatgcctaatcatcagaggctcgcctggtactgggaagtccttggcaacagggattattgctagagccatagcagacaagtaccactccagtgtgtattccttacctccagacccagaccactttgacggatacaaacaacagatcgtcactgttatggacgacctatgccaaaacccagacgggaaagacatgtcactattttgtcagatggtctccacagtggattttataccgcctatggcatctctggaggagaagggagtctcattcacctccaagtttgtgattgcctccactaacgccagtaacatcatagtgccaacagtctcggattcagatgccattcgtcgccggttctttatggactgcgatattgaggtgaccgattcctataagacagagctgggcagacttgatgcagggagagcagccaggctgtgctctgagaacaacactgcaaactttaaacggtgcagtccattagtctgtgggaaagcaatccagcttagggataggaagtccaaggtgagatacagtgtggacacggtagtgagtgaacttatcagggagtataacaacagatcagttattgggaacaccattgaagctcttttccaaggaccccctaaatttagaccaataaggattagcttagaggagaagcccgcacctgatgctattagtgacttattagctagtgttgatagtgaagaggttcgccaatactgtagagatcagggatggattgtacctgattctcccaccaacgttgagcgccacttgaatagagctgtcttgattatgcagtctgtagccaccgtggtagcagttgtgtcccttgtttacgtcatctacaagttgttcgccggttttcaaggagcatattccggcgcccccaagcaaacactcaagaaaccagtgctgcgcacggcaactgtgcaggggccgagcttggacttcgccctatctctacttaggaggaacattaggcaggtccaaaccgaccagggccactttacaatgttaggagtgcgagaccgcttggctgtgctccccagacactcccaaccaggaaagaccatctgggttgaacacaaattagtgaagatcgtagatgctgtggagttagtagacgaacaaggggttaacttagagctcacactggtaacgcttgatactaacgaaaaatttagagacatcacaagattcataccagaaacaattagtcctgctagtgatgccactttagttataaatactgaacatatgcccagtatgtttgtgccagttggagatgtggtccagtatgggtttttgaaccttagtggtaagcccactcacaggactatgatgtacaatttcccaacaaaagcaggacagtgtggtggtgttgtgactgccgtgggtaaagtgattgggatccacattggtggcaacggtaggcaaggtttctgcgctgccctgaagaggggatacttttgcagtgaacaaggtgagatccaatggatgaagcccaacaaagaaactggcaggttgaacatcaacggacctactcgcactaagcttgaaccaagtgtctttcacgatgtgttcgaaggcactaaagagccagcagtgctgactagtaaagacccaaggctggaagttgactttgaacaggctcttttttcaaaatacgtggggaacacgcttcatgaacccgacgagtttgtcaaggaggcggccttacattatgccaaccaactcaagcagttagatatcaagaccaccaagatgagcatggaggatgcatgttacggcacagagaacctggaagctatagatcttcacacaagtgcaggatatccatacagtgcactaggcatcaagaaaaaggacattttggatccaacaactcgcgatgtcagcaagatgaaattctacatggacaagtatgggttggatctaccgtactctacttatgttaaagatgaacttagggccatcgacaagatcaagaaagggaagtctcgtctcatagaagcgagcagtctaaatgactcagtgtacttgagaatgacatttgggcacctttatgaagctttccacgccaacccaggtacaatcactggttcagctgttgggtgtaacccagatgtgttctggagcaagttaccaattctacttccaggatcgcttttcgcgtttgactactcggggtatgacgctagtctcagcccagtgtggttcagggcgctggagatagtcctgcgggaaattggatactccgaagacgcagtgtctctcatagaagggatcaatcacacccatcatgtgtaccgcaataaaacttattgtgttcttgggggaatgccctcaggttgctcaggcacctccattttcaactcgatgatcaacaatatcattattagaacactcctgattaaaacattcaaagggatagatctagatgaactgaacatggtggcctacggggatgatgtgttggctagttaccccttcccaattgactgtctggagttggcaagaacaggcaaggagtatggtctaactatgacccctgccgacaagtcaccctgctttaatgaggttacatgggagaatgccactttcttgaagagaggattcttgcctgatcatcaattcccgtttctcatccaccctacgatgccaatgagggagattcacgaatccattcgttggaccaaagatgcacgaagtactcaagatcacgtgcgctccctctgcttattagcatggcacaacgggaaagaggagtatgaaaaatttgtgagtgcaatcagatcagttccaattggaaaagcattggctataccaaattatgagaatctgagaagaaattggctcgaattgttttaaatttacagtttgtaactgaaccccaccagtaatctggtcgcgttaatgactggtgggggtaaatttgttataaccagaatagc"

#e.g. python DelMapper.py /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0087_Miseq_Capsid_Deletion_Library_EV71_WalkerOrr/QVEU_Seq0087_WO_WB_04_Capsiddel_Replicationdel_DIMPLE/230820_M02211_0114_000000000-KD2BY/demux/output_test.sam 9 /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0087_Miseq_Capsid_Deletion_Library_EV71_WalkerOrr/QVEU_Seq0087_WO_WB_04_Capsiddel_Replicationdel_DIMPLE/230820_M02211_0114_000000000-KD2BY/demux/nsPs_9D_PTD.csv

#Imports 

import pandas as pd
import re
import numpy as np
import sys
from Bio import Seq
import difflib
#import pysam

#Functions
#  posParse:
#    Function to parse the cigar list in forward or reverse orientation based on alignment. Finds 'opening' (5') position of deletion on template. 
#
def posParse(cigarList):
	offset=0
	delList=[]
	for i in cigarList:
		if i[-1]=="M":#if match
			offset+=int(i[:-1])#extract length of cigar matches
		elif i[-1]=="I":#if insertion
			pass#no offset
		elif i==str(delSize)+"D":
			delList.append(offset)#if matches deletion of desired size
			offset+=int(delSize) # change PD 
		elif i[-1]=="D":
			offset+=int(i[:-1])
	return(delList)

#Functions
#  posParse:
#    Function to parse the cigar list in forward or reverse orientation based on alignment. Finds 'opening' (5') position of deletion on template. 
#    
def parselist(difflist):
	nlist=[]
	n=-1
	for i in difflist:
		if i.startswith(" "):
			n+=1
		elif i.startswith("+"):
			n=n
		elif i.startswith("-"):
			n+=1
		nlist=nlist+[n]
	return(nlist)
 #Functions
#  alignstrings():
#    Function to parse the cigar list in forward or reverse orientation based on alignment. Finds 'opening' (5') position of deletion on template. 
#   
def alignstrings(ref,query):
	d = difflib.Differ()
	aligned=list(d.compare(ref,query))
	#print(aligned)
	nlist=parselist(aligned)
	#print(nlist)
	hits=np.argwhere([i.startswith("-") for i in aligned])
	#print(hits)
	return([nlist[i[0]] for i in hits])        

#MAIN
if __name__=="__main__":
	#Arguments
	infile=sys.argv[1] #read infile argument
	delSize=sys.argv[2] #size of deletions (nt)
	outfile=sys.argv[3] #outfile (.csv)
	
	#read lines in input sam files
	with open(infile, "r") as IF:
		with open(outfile,'w') as OF:
			#linecount=0
			for line in IF:
				#linecount+=1
				if line[0]!="@":
					#print("".join(["*"]*int(linecount/100000)),end="\r")
					currentL=line.strip().split("\t")
					cigar_list=re.findall("[0-9]*[MDIX]", currentL[5])
					anchorPos=int(currentL[3])+1
					reference = Seq.translate(TEMPLATE[(anchorPos-int(currentL[3])%3):(anchorPos-int(currentL[3])%3)+(len(currentL[9])+3)])
					peptide= Seq.translate(currentL[9][2-int(currentL[3])%3::])
					diffs= alignstrings(reference[0:len(peptide)],peptide)
					deloffset=posParse(cigar_list)
					if(deloffset) and len(diffs)>0 and len(diffs)<(1+int(int(delSize)/3)):
						#print(deloffset[0])
						hit=str(int(currentL[3])+deloffset[0])
						#print(diffs)
						delPos=[d*3+int(currentL[3])-int(currentL[3])%3 if isinstance(d,int) else 0 for d in diffs]
						currentL=currentL[:11]+[hit]+[reference]+[peptide]+[str(delPos)]+[str(delPos[0])]
						#print(currentL)
						OF.write("\t".join(currentL)+"\n")
	
	print("\nCompleted.\nOutfile ("+outfile+") written.")
