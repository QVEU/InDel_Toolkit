#!/usr/bin/env python
# coding: utf-8
'''
smelt.py

Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit

#usage: python smelt.py <inputfile (str)> <delSize (int)> <outputfile (str)> 

Designed to map deletions of defined sizes on a template sequence (plasmid or transcript). Currently hard coded to map to EV-A71 reference coding sequence. 

v0.2 version notes:
    Long inserts created an issue with mapping, and highlight a issue specific to nanopore, that a good match for the insert may be flanked by bad sequence or adapters.
    Fixed by limiting the flanking sequence to 25 for all hits (regardless of the size of the insert)

GitHub push (8/14/24):
- cleaned up comments and checked out some repairs suggested by reviewer.
- Generated readmes and example data to accompany github repo. 

'''
 
#EV-A71 template
TEMPLATE="ttaaaacagcctgtgggttgcacccacccacagggcccactgggcgctagcactctggtactgaggtacctttgtgcgcctgtttttactccccttcccccgaagtaacttagaagctgtaaatcaacgatcaatagcaggtgtggcacaccagtcataccttgatcaagcacttctgtttccccggactgagtatcaataggctgctcgcgcggctgaaggagaaaacgttcgttacccgaccaactacttcgagaagcttagtaccaccatgaacgaggcagggtgtttcgctcagcacaaccccagtgtagatcaggctgatgagtcactgcaacccccatgggcgaccatggcagtggctgcgttggcggcctgcccatggagaaatccatgggacgctctaattctgacatggtgtgaagagcctattgagctagctggtagtcctccggcccctgaatgcggctaatcctaactgcggagcacatgctcacaaaccagtgggtggtgtgtcgtaacgggcaactctgcagcggaaccgactactttgggtgtccgtgtttccttttattcctatattggctgcttatggtgacaatcaaagagttgttaccatatagctattggattggccatccggtgtgcaacagggcaattgtttacctatttattggttttgtaccattatcactgaagtctgtgatcactctcaaattcattttgaccctcaacacaatcaaacatgggctcacaagtgtccacacaacgctccggttcacacgaaaactctaactcagctaccgagggttccactataaactatactaccattaattactataaagattcctatgccgccacagcaggtaagcagagccttaagcaggacccagacaagtttgcaaatcctgtcaaagacatcttcactgaaatggcagcgccattaaaatctccatctgctgaggcatgtggttacagcgatcgggtggcacaattaactattggcaattctaccatcactacgcaagaagcagcaaacatcatagttggctatggtgagtggccttcctactgttcggactctgatgctactgcagtggacaaaccaacgcgcccagatgtttcggtgaataggttttacacattggacacaaaattgtgggagaaatcatccaaggggtggtactggaaattcccggatgtgttaactgaaaccggggtctttggtcaaaatgcacagttccactacctctatcggtcagggttctgcattcacgtgcagtgcaatgctagtaagttccaccaaggagcactcctagtcgctgtcctcccagagtatgtcattgggacagtggcaggtggcacagggacggaggatagccaccccccttataagcagactcaacccggtgctgatggcttcgaattgcaacacccgtacgtgcttgatgctggcattccaatatcacaattaacagtgtgcccacatcagtggattaatttgaggaccaacaattgtgccacaataatagtgccgtacataaacgcactaccctttgattctgccttgaaccattgtaactttggtctgctggttgtgcctattagcccgttagattatgaccaaggtgcgacgccagtgatccccattactatcactttggccccaatgtgttctgaatttgcaggccttagacaagcagttacgcaagggtttcctactgagctgaaacctggcacaaaccaatttttaaccactgacgatggcgtctcagcacccattctgccaaactttcaccccaccccgtgtatccatatacccggtgaagttagaaacttgctagagctatgccaggtggagaccattttagaggtcaacaatgtacctacgaatgccactagcttaatggagagactgcgcttcccggtctcagctcaagccgggaaaggtgagctatgtgcagtgttcagagctgaccctggacgaagtgggccatggcagtccaccttgttgggccagttgtgcgggtactacacccaatggtcaggatcactggaagtcaccttcatgttcaccgggtcctttatggctaccggcaagatgctcatagcatacacaccaccaggaggccccttacccaaggaccgggcgaccgccatgttgggcacgcacgtcatctgggactttgggctgcaatcgtctgtcacccttgtaataccatggatcagcaacactcattacagagcgcacgctcgagatggtgtgttcgactactacactacaggtttggttagcatatggtaccagacgaattatgtggttccaattggggcacccaatacagcctatataatagcattggcggcagcccagaagaacttcaccatgaagttgtgtaaggatgctagtgatatcctacagacaggcactatccagggagatagggtggcagatgtgattgagagttctataggggacagtgtgagcagagccctcacccgagctctaccggcacctaccggccaagacacacaggtaagcagccaccgattagatactggtaaagttccagcactccaagccgctgaaattggagcatcatcaaatgctagtgatgagagtatgattgagacacggtgtgttcttaattcacatagtacagctgagaccactcttgatagcttcttcagcagagcaggattagttggagagatagacctccctcttgaaggcacaaccaacccgaatgggtacgcaaactgggacatagacataacaggttacgcgcaaatgcgtagaaaggtggagctgttcacctacatgcgttttgacgcagagttcacctttgttgcatgcacccctaccgggcaagttgtcccgcaattgctccaatacatgtttgtaccacccggagcccccaagccagactccagagaatctctcgcatggcaaactgccactaatccctcagtttttgtgaagctgtcagaccccccagcacaggtttctgttccattcatgtcacctgcgagcgcctatcaatggttttatgacgggtatcccacattcggtgaacacaaacaggagaaagaccttgaatacggggcatgcccaaacaacatgatgggtacgttctcagtgcggactgtaggcacctcgaagtccaagtacccattggtgatcaggatttacatgaggatgaagcacgtcagggcgtggatacctcgcccaatgcgtaaccagaactatctattcaaagccaacccaaattatgctggtaattttattaaaccaactggtgccagtcgcacagcaatcaccaccctcgggaaatttggacagcagtccggagctatctacgtgggcaactttagagtggttaaccgccatcttgctactcataatgactgggcaaaccttgtttgggaagacagctcccgcgacttgctcgtatcatctaccactgctcaaggttgtgacacgattgctcgttgcaattgccagacaggagtgtattattgtaactcaatgagaaaacactatccggtcagtttctcgaaacccagtttgatcttcgtggaggccagcgagtattatccagctagataccagtcacatctcatgcttgcagtgggtcattcggaaccaggggattgcggtggcattcttagatgccaacatggcgtcgtagggatagtttccaccgggggaaacggcctggtggggttcgccgatgtgagggatcttctgtggttggatgatgaagccatggagcagggcgtgtctgattacattaaagggcttggagatgcttttggcatggggtttacagacgcagtgtcaagagaagttgaagcactgaaaagtcacttgatcggctcagagggtgccgtggagaagattctaaagaacttagttaaactcatctctgcgctcgtcatcgtcatcaggagtgattatgacatggtcacattgacggcaacacttgccctgatcgggtgccacgggagcccttgggcctgggttaagtcgaagacagcatcaatcttgggcataccgatggctcagaagcagagtgcctcttggttaaagaagttcaacgatgcggcgagtgccgcgaaggggcttgagtggatctccaacaaaatcagtaaatttatcgattggctcaaggagaaaatcataccggctgctaaagagaaagtcgagtttctaaacaatctaaagcaactccccttattggagaaccaaatttctaatctcgaacagtcagcagcttcgcaggaggaccttgaggcgatgtttggcaacgtgtcttatctggcccacttctgccgcaaattccaacccctctatgccacggaagcaaagagggtgtacgccctagaaaagagaatgaataattacatgcagttcaagagcaaacaccgtattgaacctgtatgcctaatcatcagaggctcgcctggtactgggaagtccttggcaacagggattattgctagagccatagcagacaagtaccactccagtgtgtattccttacctccagacccagaccactttgacggatacaaacaacagatcgtcactgttatggacgacctatgccaaaacccagacgggaaagacatgtcactattttgtcagatggtctccacagtggattttataccgcctatggcatctctggaggagaagggagtctcattcacctccaagtttgtgattgcctccactaacgccagtaacatcatagtgccaacagtctcggattcagatgccattcgtcgccggttctttatggactgcgatattgaggtgaccgattcctataagacagagctgggcagacttgatgcagggagagcagccaggctgtgctctgagaacaacactgcaaactttaaacggtgcagtccattagtctgtgggaaagcaatccagcttagggataggaagtccaaggtgagatacagtgtggacacggtagtgagtgaacttatcagggagtataacaacagatcagttattgggaacaccattgaagctcttttccaaggaccccctaaatttagaccaataaggattagcttagaggagaagcccgcacctgatgctattagtgacttattagctagtgttgatagtgaagaggttcgccaatactgtagagatcagggatggattgtacctgattctcccaccaacgttgagcgccacttgaatagagctgtcttgattatgcagtctgtagccaccgtggtagcagttgtgtcccttgtttacgtcatctacaagttgttcgccggttttcaaggagcatattccggcgcccccaagcaaacactcaagaaaccagtgctgcgcacggcaactgtgcaggggccgagcttggacttcgccctatctctacttaggaggaacattaggcaggtccaaaccgaccagggccactttacaatgttaggagtgcgagaccgcttggctgtgctccccagacactcccaaccaggaaagaccatctgggttgaacacaaattagtgaagatcgtagatgctgtggagttagtagacgaacaaggggttaacttagagctcacactggtaacgcttgatactaacgaaaaatttagagacatcacaagattcataccagaaacaattagtcctgctagtgatgccactttagttataaatactgaacatatgcccagtatgtttgtgccagttggagatgtggtccagtatgggtttttgaaccttagtggtaagcccactcacaggactatgatgtacaatttcccaacaaaagcaggacagtgtggtggtgttgtgactgccgtgggtaaagtgattgggatccacattggtggcaacggtaggcaaggtttctgcgctgccctgaagaggggatacttttgcagtgaacaaggtgagatccaatggatgaagcccaacaaagaaactggcaggttgaacatcaacggacctactcgcactaagcttgaaccaagtgtctttcacgatgtgttcgaaggcactaaagagccagcagtgctgactagtaaagacccaaggctggaagttgactttgaacaggctcttttttcaaaatacgtggggaacacgcttcatgaacccgacgagtttgtcaaggaggcggccttacattatgccaaccaactcaagcagttagatatcaagaccaccaagatgagcatggaggatgcatgttacggcacagagaacctggaagctatagatcttcacacaagtgcaggatatccatacagtgcactaggcatcaagaaaaaggacattttggatccaacaactcgcgatgtcagcaagatgaaattctacatggacaagtatgggttggatctaccgtactctacttatgttaaagatgaacttagggccatcgacaagatcaagaaagggaagtctcgtctcatagaagcgagcagtctaaatgactcagtgtacttgagaatgacatttgggcacctttatgaagctttccacgccaacccaggtacaatcactggttcagctgttgggtgtaacccagatgtgttctggagcaagttaccaattctacttccaggatcgcttttcgcgtttgactactcggggtatgacgctagtctcagcccagtgtggttcagggcgctggagatagtcctgcgggaaattggatactccgaagacgcagtgtctctcatagaagggatcaatcacacccatcatgtgtaccgcaataaaacttattgtgttcttgggggaatgccctcaggttgctcaggcacctccattttcaactcgatgatcaacaatatcattattagaacactcctgattaaaacattcaaagggatagatctagatgaactgaacatggtggcctacggggatgatgtgttggctagttaccccttcccaattgactgtctggagttggcaagaacaggcaaggagtatggtctaactatgacccctgccgacaagtcaccctgctttaatgaggttacatgggagaatgccactttcttgaagagaggattcttgcctgatcatcaattcccgtttctcatccaccctacgatgccaatgagggagattcacgaatccattcgttggaccaaagatgcacgaagtactcaagatcacgtgcgctccctctgcttattagcatggcacaacgggaaagaggagtatgaaaaatttgtgagtgcaatcagatcagttccaattggaaaagcattggctataccaaattatgagaatctgagaagaaattggctcgaattgttttaaatttacagtttgtaactgaaccccaccagtaatctggtcgcgttaatgactggtgggggtaaatttgttataaccagaatagc"


#Imports 

import pandas as pd
import re
import numpy as np
import sys
from Bio import Seq
import difflib

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
#  parselist: takes list of output of difflib function and tabulates the position by counting the insertions and deletions relative to the reference. 

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
			#For each line in SAM
			for line in IF:
				#skip comments and...
				if line[0]!="@":
					#split SAM by tabs
					currentL=line.strip().split("\t")
					# Use REGEX to find all chunks of cigar string (element 5 in SAM)
					cigar_list=re.findall("[0-9]*[MDIX]", currentL[5])
					
					# Adjust for 0 index
					anchorPos=int(currentL[3])+1
					
					# translate reference
					reference = Seq.translate(TEMPLATE[(anchorPos-int(currentL[3])%3):(anchorPos-int(currentL[3])%3)+(len(currentL[9])+3)])
					
					# translate read adjusting for reading frame
					peptide= Seq.translate(currentL[9][2-int(currentL[3])%3::])
					
					# find differences between strings using alignstrings()
					diffs= alignstrings(reference[0:len(peptide)],peptide)
					
					#use posParse() to 
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
