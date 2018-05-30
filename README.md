# sampleBiMarkers
Code for sampling bi-allelic markers from samples or populations for running PhyloNet's MLE_biMarkers. Feel free to use this code however you see fit.

Requires Python 3. 

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
		-a,--allowN		: Toggle on to allow N to be sampled
		-g,--allowG		: Toggle on to allow gap characters to be sampled
		-G,--keepG		: Toggle on to NOT treat gap characters as missing data for -N, -n options
		-m,--allowM		: Toggle on to allow loci that are monomorphic as a consequence of random sampling
			-When sampling a small number of alleles, it is possible to lose all variation.
			-Use this option to turn off the monomorphic filter that is applied AFTER sampling alleles
		-h,--help	: Displays help menu

		Note that this script will not sample Ns or gap characters by default. Both are treated as missing data. Filter accordingly.
```

