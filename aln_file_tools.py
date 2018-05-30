#!/usr/bin/python

import re
import sys
import os

"""Functions for parsing and manipulating sequence alignment files
Functions by Zach Zbinden and Tyler Chafin"""


#Function to write an alignment as DICT to NEXUS
def dict2nexus(nex, aln):
	with open(nex, 'w') as fh:
		try:
			slen = getSeqLen(aln)
			header = "#NEXUS\n\nBegin data;\nDimensions ntax=" + str(len(aln)) + " nchar=" + str(slen) + ";\n"
			header = header + "Format datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n\n"
			fh.write(header)
			for seq in aln:
				sline = str(seq) + " " + aln[seq] + "\n"
				fh.write(sline)
			last = ";\nEnd;\n"
			fh.write(last)
		except IOError:
			print("Could not read file ",nex)
			sys.exit(1)
		finally:
			fh.close()

#Reads an alignment as FASTA. returns dict of lists
#key = sample name (FASTA header)
#list = split sequence by character
#handles interleaved or non-interleaved FASTA files
def readFastaAlign(fas):
	if not fileCheck(fas):
		raise FileNotFoundError("Fatal exception, file %s not found."%fas)

	fh = open(fas)
	try:
		ret = dict()
		with fh as file_object:
			contig = ""
			seq = ""
			for line in file_object:
				line = line.strip()
				if not line:
					continue
				line = line.replace(" ","")
				#print(line)
				if line[0] == ">": #Found a header line
					#If we already loaded a contig, yield that contig and
					#start loading a new one
					if contig:
						ret[contig] = seq.split()
						contig = "" #reset contig and seq
						seq = ""
					contig = (line.replace(">",""))
				else:
					seq += line
		#Iyield last sequence, if it has both a header and sequence
		if contig and seq:
			ret[contig] = seq.split()
		return(ret)
	finally:
		fh.close()

#Function to read a phylip file. Returns dict (key=sample) of lists (sequences divided by site)
def readPhylip(phy):
	if os.path.exists(phy):
		with open(phy, 'r') as fh:
			try:
				num=0
				ret = dict()
				for line in fh:
					line = line.strip()
					if not line:
						continue
					num += 1
					if num == 1:
						continue
					arr = line.split()
					ret[arr[0]] = list(arr[1])
				return(ret)
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)


#Write FASTA from pandas df where col1 is index, col2 is sequence
#seqs must be a pandas df
def writeFasta(seqs, fas):
	file_object = open(fas, "w")
	#Write seqs to FASTA first
	#Assumes that a[0] is index, a[1] is id, and a[2] is sequence
	for a in seqs.itertuples():
		name = ">id_" + str(a[1]) + "\n"
		seq = a[2] + "\n"
		file_object.write(name)
		file_object.write(seq)
	file_object.close()
