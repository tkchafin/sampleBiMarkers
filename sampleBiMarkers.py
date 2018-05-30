#!/usr/bin/python

import re
import sys
import os
import getopt
import operator
import random
import aln_file_tools as aln
import seq_tools as seq
import collections

"""NOTE: When writing this, I prioritized speed of writing the code (e.g. minimizing time
to getting the input file I wanted), and so it is very inefficient and does not
implement robust error handling. Sorry. -TKC"""

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

		#Remove samples from pop_assign that do not have data
		pop_assign = aln.cleanPopmap(pop_assign, seqs.keys())

		#Make dict of dicts that splits by population, only retaining pops/samples from the popmap.
		#Get unique pop names using one of the worst python lines ever written
		pops = dict()
		for k in set(pop_assign.values()):
			pops[k] = dict()

		#Remove pops listed as excluded
		if params.exclude:
			print("Excluding populations:", ", ".join(params.exclude))
			for exc in params.exclude:
				if exc in pops:
					del pops[exc]
		if params.include:
			print("Only keeping populations:", ", ".join(params.include))
			for pop in list(pops):
				if pop not in params.include:
					del pops[pop]

		#make sure we didn't throw out all populations...
		if len(list(pops)) < 1:
			print("Oops! No populations remaining. Check that popmap sample names match those in your data file, or that selections using --include or --exclude are correct! :)")
			sys.exit(1)

		#Capture samples for each pop in a dict of dicts
		#Note that this removes samples for which we have data but no pop assignment
		alen = seq.getSeqLen(seqs)
		for assigned in pop_assign:
			if pop_assign[assigned] in pops:
				pops[pop_assign[assigned]][assigned] = seqs[assigned]
		seqs.clear()
		#debugging print
		# for pop in pops.keys():
		# 	print("\n",pop)
		# 	for ind in pops[pop]:
		# 		print(ind, len(pops[pop][ind]), "nucleotides")

		#Make 2D list to remove columns failing the globalN filter
		bad_columns = list() #list of column numbers to delete

		#For each pop dict, make 2D list to remove columns failing popN filter
		print("Found",alen,"nucleotide columns in the dataset!")
		columns = [[]for i in range(alen)] #2D array of global data
		for pop, data in pops.items():
			for sample, sequence in data.items():
				for i, nuc in enumerate(sequence):
					columns[i].append(nuc)

		#Remove columns with >globalN ambiguous bases, or those which are non-biallelic
		print("Removing non-biallelic columns...")
		if not params.keepG:
			print("Checking N content (counting gaps as missing data)...")
		else:
			print("Checking N content (NOT counting gaps as missing data)...")
		for i, col in enumerate(columns):
			if not seq.isBiallelic(col):
				bad_columns.append(i)
			if not params.keepG:
				if seq.checkNGcontent(col, params.globalN): #counts gaps as Ns
					bad_columns.append(i)
			else:
				if seq.checkNcontent(col, params.globalN): #counts gaps as Ns
					bad_columns.append(i)
		bad_columns=sorted(set(bad_columns))
		columns.clear()

		#Now get bad columns WITHIN pops for localN failing sites
			#New 2D structure: dict of pops, each pop a list of lists.
			#Lists in each pop are COLUMNS of nucleotide data
			#From here on, we don't track which alleles belong to which individual
			#the end product will be alleles sampled at random from each population
		popFinal = dict()
		for pop, data in pops.items():
			columns = [[]for i in range(alen)]
			for sample, sequence in data.items():
				for i, nuc in enumerate(sequence):
					columns[i].append(nuc)
			for i, col in enumerate(columns):
				if not params.keepG:
					if seq.checkNGcontent(col, params.globalN): #counts gaps as Ns
						bad_columns.append(i)
				else:
					if seq.checkNcontent(col, params.globalN): #counts gaps as Ns
						bad_columns.append(i)
			popFinal[pop] = columns
		#Re-sort and unique the bad_columns list
		bad_columns=sorted(set(bad_columns))
		print("Deleting",len(bad_columns),"columns!")
		#For each pop, remove blacklisted columns
		for pop in popFinal:
			for col in sorted(bad_columns,reverse=True): #reverse sorted so subsequent deletes aren't thrown off
				#print(col,":",columns[col])
				del popFinal[pop][col]

		#Finally, sample alleles from pops.
		outputFinal = collections.OrderedDict()
		outputAssign = collections.OrderedDict()
		for pop in popFinal:
			for i in range(0,params.sample):
				name=pop + "_" + str(i)
				outputAssign[name] = pop
				outputFinal[name] = list()
			#print(outputAssign)
			for nucs in popFinal[pop]:
				#sample alleles at random, without replacement
				sampled = seq.sampleAlleles(nucs, params.sample, params.allowN, params.allowG)
				#check that we sampled enough alleles
				if len(sampled) < params.sample :
					print("Uh oh! Not enough valid alleles in pop",pop,"! Unfortunately, my developer was too lazy to write better error checking! More stringeng --maxN and --popN filtering, or sample less alleles!!")
					sys.exit(1)

				#insert alleles into outputFinal
				for i in range(0,params.sample):
					name=pop + "_" + str(i)
					outputFinal[name].append(sampled[i])
		popFinal.clear()

		#Convert to numeric format...
		outputFinal_formatted = collections.OrderedDict()
		for ind in outputFinal:
			#print(ind,", Nucs:",outputFinal[ind][0:10])
			outputFinal_formatted[ind] = list()

		if not params.allowM:
			print("Now applying a final check for monomorphic loci...")
		print("Converting to numeric format...")
		countM = 0
		for i in range(0,seq.getSeqLen(outputFinal)):
			this_column = list()
			for ind in outputFinal:
				this_column.append(outputFinal[ind][i])
			if seq.isMonomorphic(this_column):
				if not params.allowM:
					continue
			numeric = seq.nucs2numeric(this_column)
			k=0
			for ind in outputFinal:
				outputFinal_formatted[ind].append(numeric[k])
				k+=1
		outputFinal.clear()

		#output alignment to NEXUS
		#for ind in outputFinal_formatted:
			#print(ind,", Nucs:",outputFinal_formatted[ind][0:10])
		aln.dict2nexus(params.out, outputFinal_formatted)

		#Append MLE_BiMarkers skeleton command to NEXUS file
		print("Writing PHYLONET block to Nexus file...")
		with open(params.out, 'a') as fh:
			try:
				line1 = "BEGIN PHYLONET;\n"
				fh.write(line1)
				taxalist = ",".join(list(outputAssign))
				taxamap = ""
				popGroup = dict()
				for sampleID, popID in outputAssign.items():
					if popID not in popGroup:
						popGroup[popID] = list()
					popGroup[popID].append(sampleID)
				first = True
				for popID in popGroup.keys():
					if first:
						taxamap = popID + ":" + ",".join(popGroup[popID])
						first = False
					else:
						taxamap += ";" + popID + ":" + ",".join(popGroup[popID])
				#print(taxamap)
				line2 = "MLE_BiMarkers -taxa (" + taxalist + ") -tm <" + taxamap + ">;\nEND;\n"
				fh.write(line2)
			except IOError:
				print("Could not read file ",params.out)
				sys.exit(1)
			finally:
				fh.close()


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
			options, remainder = getopt.getopt(sys.argv[1:], 'f:i:ho:dp:s:N:n:x:I:agGm', \
			["input=","phylip=","phy=","out=","nohet","fasta=","popmap=","maxN=",
			"popN=","exclude=","include=", "allowG", "allowN", "keepG","allowM"])
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
		self.include = list()

		#booleans
		self.allowN=False
		self.allowG=False
		self.allowM=False
		self.keepG=False
		self.noHet=False

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
			elif opt in ('I','include'):
				self.include = arg.split(",")
			elif opt in ("a", "allowN"):
				self.allowN=True
			elif opt in ("g", "allowG"):
				self.allowG=True
			elif opt in ("G", "keepG"):
				self.keepG=True
			elif opt in ("m", "allowM"):
				self.allowM=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.phylip and not self.fasta:
			self.display_help("Error: Missing required alignment file (--fasta or --input)")
		if not self.popmap:
			self.display_help("Error: Missing required popmap file (-p, --popmap)")
		if self.include and self.exclude:
			self.display_help("Don't use both --include and --exclude.")
		if self.sample < 1:
			self.display_help("Error: --sample must be greater than 1.")


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
		-x,--exclude	: List of pops to exclude (format: -x "Pop1,Pop2,Sample4...")
		-I,--include	: List of pops to include (removing all others)
		-a,--allowN		: Toggle on to allow N to be sampled
		-g,--allowG		: Toggle on to allow gap characters to be sampled
		-G,--keepG		: Toggle on to NOT treat gap characters as missing data for -N, -n options
		-m,--allowM		: Toggle on to allow loci that are monomorphic as a consequence of random sampling
			-When sampling a small number of alleles, it is possible to lose all variation.
			-Use this option to turn off the monomorphic filter that is applied AFTER sampling alleles
		-h,--help	: Displays help menu

		Note that this script will not sample Ns or gap characters by default. Both are treated as missing data. Filter accordingly.

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
