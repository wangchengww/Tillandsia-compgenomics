#!/usr/bin/env python
import sys
import collections
import argparse
import copy
import math
import random
import re


class FastaReader:
	"""
	A light-weight fasta reader;
	returns a tuple (header, sequence)
	
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
		self.__prevheader=None

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()
	
	def next(self):
		line=""
		header=self.__prevheader
		seq=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":					# file is empty
				if(header is not None):
					self.__prevheader=None		# return last sequence
					return (header,seq)
				else:
					raise StopIteration		# file empty and no last sequence -> STOP
			line=line.rstrip("\n")				# somethin is in the file
			if(line.startswith(">")):			# we have a header
				line=line.lstrip(">")
				if(header is None):			# if it is the first header just set the name of the header
					header=line
				else:
					self.__prevheader=line	# if it is any other header, set the previous to the current and return the sequence
					return(header,seq)
			else:
				seq+=line				# normal line, add to sequence
	
	
	@classmethod
	def readAllTuples(cls,file):
		ft=[]
		for n,s in FastaReader(file):
			ft.append((n,s))
		return ft

class FastaWriter:
	"""
	Write the content to a fasta file
	"""
	def __init__(self,file,seqleng):
		self.__filename=file
		self.__filehandle=open(file,"w")
		self.__seqleng=seqleng
		
	def write(self,n,s):
		sl=self.__seqleng
		fh=self.__filehandle
		fh.write(">"+n+"\n")
		c=0
		while(c<len(s)):
			fh.write(s[c:c+sl]+"\n")
			c+=sl

	def close(self):
		self.__filehandle.close()
	
	@classmethod
	def write_all(cls,file,fastaentries):
		fw=FastaWriter(outputFile)
		for n,s in fastaentries:
			fw.write(n,s)
		fw.close()
		return 1


def read_breakpoints(bpfile):
    bph=collections.defaultdict(lambda:[])
    for l in open(bpfile):
        a=l.rstrip("\n").split(":")
        chrm,pos=a[0],int(a[1])
        bph[chrm].append(pos)
    
    # sort descending largest first
    for chrm in bph.keys():
        vals=bph[chrm]
        vals=sorted(vals,key=lambda i:-i)
        bph[chrm]=vals
    return bph




parser = argparse.ArgumentParser(description="""           
Description
-----------
    introduce breaks into fasta""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 2.7


""")


parser.add_argument("--input", type=str, required=True, dest="inp", default=None, help="A fasta file, where breakpoints should be introduced ")
parser.add_argument("--output", type=str, required=True, dest="out", default=None, help="The output file; a fasta file")
parser.add_argument("--breakpoints", type=str, required=True, dest="breaks", default=None, help="a file with breakpoints, each line a bp with format 'chr:pos' ")
args = parser.parse_args()

bphash=read_breakpoints(args.breaks)
fr=FastaReader(args.inp)
fw=FastaWriter(args.out,60)

for header,seq in fr:
    if header not in bphash:
        fw.write(header,seq)
        continue
    bpsForChr=bphash[header]
    
    workseq=seq
    
    for bp in bpsForChr:
        abspatz=workseq[bp:]
        workseq=workseq[:bp]
        fw.write(header+str(bp),abspatz)
    fw.write(header+"rest",workseq)

fr.close()
fw.close()
    


