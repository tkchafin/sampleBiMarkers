#!/usr/bin/python

import re
import sys
import os

"""Tools for manipulating sequence data by Tyler Chafin"""

#Function takes biallelic list of nucleotides and converts to numeric
#0 = major allele
#1 = minor allele
#2 = het
#? = - or N
def nucs2numeric(nucs):
	if isBiallelic(nucs):
		#print(nucs)
		ret = list()
		counts = {"A":0, "G":0, "C":0, "T":0}
		#find major allele
		for nuc in nucs:
			if nuc not in ("-", "N"):
				for exp in get_iupac_caseless(nuc):
					counts[exp] += 1
		#sort dict, to list of tuples (b/c dicts are orderless, can't keep as dict)
		sorted_x = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
		majA = sorted_x[0][0]
		minA = sorted_x[1][0]
		het = reverse_iupac(''.join(sorted(set([majA, minA])))) #get het code
		#print(majA, minA, het)

		for nuc in nucs:
			nuc = nuc.upper()
			if nuc == majA:
				ret.append("0")
			elif nuc == minA:
				ret.append("1")
			elif nuc == het:
				ret.append("2")
			elif nuc == "-":
				ret.append("-")
			else:
				ret.append("?")

		return(ret)
	else:
		print("Warning: Data is not biallelic:",nucs)
		return(None)

#Function takes biallelic list of nucleotides and converts to numeric
#0 = major allele
#1 = minor allele
#2: Randomly samples heterozygous sites as 0 or 1
def nucs2numericNohet(nucs):
	if isBiallelic(nucs):
		#print(nucs)
		ret = list()
		counts = {"A":0, "G":0, "C":0, "T":0}
		#find major allele
		for nuc in nucs:
			if nuc not in ("-", "N"):
				for exp in get_iupac_caseless(nuc):
					counts[exp] += 1
		#sort dict, to list of tuples (b/c dicts are orderless, can't keep as dict)
		sorted_x = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
		majA = sorted_x[0][0]
		minA = sorted_x[1][0]
		het = reverse_iupac(''.join(sorted(set([majA, minA])))) #get het code
		#print(majA, minA, het)

		for nuc in nucs:
			nuc = nuc.upper()
			if nuc == majA:
				ret.append("0")
			elif nuc == minA:
				ret.append("1")
			elif nuc == het:
				ret.append(random.randint(0,1))
			elif nuc == "-":
				ret.append("-")
			else:
				ret.append("?")

		return(ret)
	else:
		print("Warning: Data is not biallelic:",nucs)
		return(None)

#Function to translate a string of bases to an iupac ambiguity code
def reverse_iupac(char):
	char = char.upper()
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N'
	}
	return iupac[char]

#Function takes a list of nucleotides, and returns True if the column is biallelic
#ignores gaps and Ns
#expands uipac codes using a call to external function
def isBiallelic(nucs):
	expanded = list()
	for nuc in nucs:
		if nuc not in ("-", "N"):
			for exp in get_iupac_caseless(nuc):
				expanded.append(exp)
	uniq_sort = sorted(set(expanded))
	if len(uniq_sort) != 2:
		#print(nucs)
		#print(uniq_sort, len(uniq_sort))
		return(False)
	else:
		return(True)

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
	if char.islower():
		char = char.upper()
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["A", "C", "G", "T"],
		"-"	: ["-"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: ["C","G","T"],
		"D"	: ["A","G","T"],
		"H"	: ["A","C","T"],
		"V"	: ["A","C","G"]
	}
	ret = iupac[char]
	return ret



#Goes through a dict of sequences and get the alignment length
def getSeqLen(aln):
	length = None
	for key in aln:
		if not length:
			length = len(aln[key])
		else:
			if length != len(aln[key]):
				print("getSeqLen: Alignment contains sequences of multiple lengths.")
	return(length)
