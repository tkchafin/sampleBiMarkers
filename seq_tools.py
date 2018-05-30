#!/usr/bin/python

import re
import sys
import os
import random
import operator

"""Tools for manipulating sequence data by Tyler Chafin"""

#function to sample non-N alleles without replacement
def sampleAlleles(geno, n, allowN, allowG):
	expanded = list()
	for gen in geno:
		for nuc in get_iupac_caseless_diploid(gen):
			if not allowN and nuc in ("N", "n"):
				continue
			if not allowG and nuc == "-":
				continue
			expanded.append(nuc)
	if n > len(expanded):
		return(expanded)
	else:
		ret = list()
		for sampled in sampleList(expanded, n):
			ret.append(sampled)
		return(ret)


#utility function to randomly sample w/out replacement from a list
def sampleList(stuff, r):
	rand = random.random
	n = len(stuff)
	pop = len(stuff)
	for samp in range(r, 0, -1):
		cumprob=1.0
		x=rand()
		while x < cumprob:
			cumprob -= cumprob*samp/pop
			pop-=1
		yield stuff[n-pop-1]


#function to check if N content in list of diploid genotypes (1 element = 1 pair of alleles)
#is greater than a given threshold
#Returns TRUE if N+gap content is too high
#this version treats Ns and gaps as the same
def checkNGcontent(nucs, threshold):
	nucs = [x.upper() for x in nucs]
	twoN = len(nucs)*2
	Ncontent = ((nucs.count("N")*2) + (nucs.count("-")*2))
	#print("Ncontent is: ",float(Ncontent/twoN), ":", nucs)
	if float(Ncontent/twoN) >= float(threshold):
		return True
	else:
		return False

#function to check if N content in list of diploid genotypes (1 element = 1 pair of alleles)
#is greater than a given threshold
#Returns TRUE if N+gap content is too high
#this version only counts N content
def checkNcontent(nucs, threshold):
	nucs = [nucs.upper() for x in nucs]
	twoN = len(nucs)*2
	Ncontent = (nucs.count("N")*2)
	if float(Ncontent/twoN) >= float(threshold):
		return True
	else:
		return False

#Function takes biallelic list of nucleotides and converts to numeric
#0 = major allele
#1 = minor allele
#2 = het
#? = - or N
def nucs2numeric(nucs):
	ret = list()
	counts = {"A":0, "G":0, "C":0, "T":0}
	#find major allele
	for nuc in nucs:
		nuc = nuc.upper()
		if nuc not in ("-", "N", "n"):
			counts[nuc]+=1

	if isBiallelic(nucs):
		#print(nucs)
		#sort dict, to list of tuples (b/c dicts are orderless, can't keep as dict)
		sorted_x = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
		majA = sorted_x[0][0]
		minA = sorted_x[1][0]
		#print(majA, minA, het)

		for nuc in nucs:
			nuc = nuc.upper()
			if nuc == majA:
				ret.append("0")
			elif nuc == minA:
				ret.append("1")
			elif nuc == "-":
				ret.append("-")
			else:
				ret.append("?")

		return(ret)
	elif isMonomorphic(nucs):
		for nuc in nucs:
			ret.append("0")
		return(ret)
	else:
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

#Function takes a list of nucleotides, and returns True if the column is monomorphic
#ignores gaps and Ns
#expands uipac codes using a call to external function
def isMonomorphic(nucs):
	expanded = list()
	for nuc in nucs:
		if nuc not in ("-", "N"):
			for exp in get_iupac_caseless(nuc):
				expanded.append(exp)
	uniq_sort = sorted(set(expanded))
	if len(uniq_sort) != 1:
		#print(nucs)
		#print(uniq_sort, len(uniq_sort))
		return(False)
	else:
		return(True)


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

#Function to split character to IUPAC codes, returning diploid genotype
def get_iupac_caseless_diploid(char):
	if char.islower():
		char = char.upper()
	iupac = {
		"A"	: ["A", "A"],
		"G"	: ["G", "G"],
		"C"	: ["C", "C"],
		"T"	: ["T", "T"],
		"N"	: ["N", "N"],
		"-"	: ["-","-"],
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
