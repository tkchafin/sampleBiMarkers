# sampleBiMarkers
Code for sampling bi-allelic markers from samples or populations for running PhyloNet's MLE_biMarkers. Feel free to use this code however you see fit.

Code has been tested and seems to create valid input files. Requires Python 3. 

### How to use
To view the help menu, call the script using the <-h,--help> argument:

Which will display:

    ./sampleBiMarkers.py -h

Which will display: 
```
uaf90168:sampleBiMarkers tkchafin$ python3 sampleBiMarkers.py -h

Exiting because help menu was called.

sampleBiMarkers.py

Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu

Usage:  sampleBiMarkers.py -p /path/to/phylip 

Description: Script for sampling bi-allelic markers for running PhyloNet MLE_Bimarkers

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
		-a,--allowN	: Toggle on to allow N to be sampled
		-g,--allowG	: Toggle on to allow gap characters to be sampled
		-G,--keepG	: Toggle on to NOT treat gap characters as missing data for -N, -n options
		-m,--allowM	: Toggle on to allow loci that are monomorphic as a consequence of random sampling
			-When sampling a small number of alleles, it is possible to lose all variation.
			-Use this option to turn off the monomorphic filter that is applied AFTER sampling alleles
		-h,--help	: Displays help menu

		Note that this script will not sample Ns or gap characters by default. Both are treated as missing data. Filter accordingly.
```
These options can be used to change the behaviour of the script. By default, sampleBiMarkers.py will treat gaps and Ns as missing data, and delete columns with too high of a proportion globally across all samples threshold controlled using <-N,--maxN> or containing too much missing data within any sampled population <-n,--popN>. You can have gaps treated as valid characters by using the <-G,--kepG> flag. After filtering for missing data, the script will sample <-s,--sample> number of alleles from each population, not sampling any missing data unless <-a,--allowN> or <-g,--allowG> are used. 

The final output is a NEXUS file formatted numerically (where 0=major allele; 1=minor allele; ?=N; -=gap), with a draft command block for the PhyloNet MLE_BiMarkers command. You will need to add any additional arguments or settings to this line (e.g. if you want to toggle on the pseudolikelihood calculation, add the flag "-pseudo" after the "MLE_BiMarkers"). 

Note that this script was written quickly, and is not terribly efficient nor does it implement robust error checking/reporting. If you have any issues with it, please email me directly at tkchafin@uark.edu

### Input requirements
Alignments of SNPs can be input using PHYLIP or FASTA format. You will need to specify how these samples are partitioned into species/populations/groups/etc using a tab-delimited 'population map', which should be formatted where column1 is the sample ID (exactly matching a header from the PHYLIP or FASTA file) and column 2 is a population ID:
```
SampleA	Population1
SampleB	Population1
SampleC	Population2
...
...
...
```
You can control which populations are included by sampleBiMarkers using the --include or --exclude flags. For example, to sample 2 alleles each from only (and create output for) Population1, Population3, and Population4:

```./sampleBiMarkers.py -i example.input -s 2 --include Population1,Population3,Population4```

Be sure that your populations are spelled correctly, and separated by commas (with NO spaces). 

One helpful tip to see a short list of your population codes, is to type the following bash command, replacing $POPMAP with the name of your popmap file:

```cat $POPMAP | awk '{print $2}' | sort | uniq```
