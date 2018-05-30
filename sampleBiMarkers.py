#!/usr/bin/python

import re
import sys
import os
import getopt
import operator
import random
import aln_file_tools as aln
import seq_tools as seq

def main():
	params = parseArgs()

	seqs = dict()
	if params.phylip:
		print("Parsing phylip file...")
		#Get sequences as dict of lists
		seqs = aln.readPhylip(params.phylip)
	elif params.fasta:
		print("Parsing fasta file...")
		seqs = aln.readFastaAlign(params.fasta)
	else:
		print("No input provided.")
		sys.exit(1)

	pop_assign = dict()
	#parse popmap file for dictionary of sample assignments
	if params.popmap:
		print("Parsing popmap file...")
		pop_assign = aln.parsePopmap(params.popmap)
	else:
		print("ERROR: Popmap file must be provided.")
		sys.exit(1)

	if seqs and pop_assign:

		#Make dict of dicts that splits by population, only retaining pops/samples from the popmap.
		#Remove 'excluded' populations
		#Make 2D list to remove columns failing the globalN filter
		#For each pop dict, make 2D list to remove columns failing popN filter
		#Delete failed columns from master 2D dict
		#Finally, sample alleles from pops. 

	else:
		print("Oops! Something went wrong reading alignment or popmap file! Sorry this error message isn't more informative! >:)")
		sys.exit(1)

	# #get list of columns and list of samplenames
	# alen = getSeqLen(seqs)
	# columns = [[]for i in range(alen)]
	# names = list()
	# for key, value in seqs.items():
	# 	names.append(key)
	# 	for i, nuc in enumerate(value):
	# 		columns[i].append(nuc)
	#
	# #For each column, delete those which are not bi-allelic
	# dels=list()
	# for i, col in enumerate(columns):
	# 	if not isBiallelic(col):
	# 		dels.append(i)
	# 		#print(i,"not biallelic:",col)
	#
	# print("Deleting",len(dels),"non-biallelic columns.")
	# for col in sorted(dels,reverse=True): #reverse sorted so subsequent deletes aren't thrown off
	# 	#print(col,":",columns[col])
	# 	del columns[col]
	#
	# #Then, convert to 012 format
	# print("Converting to 012 format...")
	# formatted = [[]for i in range(alen-len(dels))]
	#
	# for i, col in enumerate(columns):
	# 	#print(col)
	# 	#print(nucs2numeric(col))
	# 	if params.nohet:
	# 		formatted[i] = nucs2numericNohet(col)
	# 	else:
	# 		formatted[i] = nucs2numeric(col)
	# 	#sys.exit()
	#
	# final_data = dict()
	# for i, samp in enumerate(names):
	# 	seqs = list()
	# 	for k,nuc in enumerate(formatted):
	# 		seqs.append(nuc[i])
	# 	final_data[samp] = "".join(seqs)
	#
	# print("Writing NEXUS output file...")
	# dict2nexus(params.out, final_data)



#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'f:i:ho:dp:s:N:n:', \
			["input=","phylip=","phy=","out=","nohet","fasta=","popmap=","maxN=","popN="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.phylip=None
		self.fasta=None
		self.popmap=None
		self.out="out.nex"
		self.nohet=False
		self.sample=1
		self.globalN=0.5
		self.popN=0.5
		self.exclude = list()

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt in ('i', 'phylip', 'input','phy'):
				self.phylip = arg
			elif opt in ('f', 'fasta'):
				self.fasta = arg
			elif opt in ('p', 'popmap'):
				self.popmap = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('s', 'sample'):
				self.sample = int(arg)
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ('d','nohet'):
				self.nohet=True
			elif opt in ('N', 'maxN'):
				self.globalN = float(arg)
			elif opt in ('n', 'popN'):
				self.popN = float(arg)
			elif opt in ('x', 'exclude'):
				self.exclude = arg.split(",")
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.phylip and not self.fasta:
			self.display_help("Error: Missing required alignment file (--fasta or --input)")
		if not self.popmap:
			self.display_help("Error: Missing required popmap file (-p, --popmap)")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nsampleBiMarkers.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-p /path/to/phylip \n")
		print ("Description: Script for sampling bi-allelic markers for running PhyloNet MLE_Bimarkers")

		print("""
	Arguments:
		INPUT FILES [REQUIRED]
		-i,--input	: Input file as PHYLIP
			-or-
		-f,--fasta	: optionally input your data as FASTA
		-p,--popmap	: Tab-delimited population map

		PARAMETERS [OPTIONAL]
		-o,--out	: Output file name <default = out.nex>
		-s,--sample	: Number of alleles to sample [default=1]
		-N,--maxN	: Maximum proportion of globally missing data allowed to drop a SNP [default=0.5]
		-n,--popN	: Maximum proportion of Ns within pop to drop SNP [default=0.5]
		-x,--exclude	: List of samples or pops to exclude (format: -x "Pop1,Pop2,Sample4...")
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
