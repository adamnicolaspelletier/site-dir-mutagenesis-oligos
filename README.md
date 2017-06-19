# site-dir-mutagenesis-oligos
Generate oligos for site directed mutagenesis with a Q5 Phusion-type home-made approach in an automated fashion

### Author: Adam-Nicolas Pelletier
### Github: https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos.git

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Information] (#information)
- [Usage](#usage)
- [In Development](#in development)
- [Known Bugs](#known bugs)
- [Support](#support)
- [Contributing](#contributing)



## Requirements

* Python 2.7+
* Python Libraries :

  1. Pandas (+ Numpy) 0.20.1+
  2. BioPython 1.69+
	

* Exon FASTA sequences for genes or transcripts of interest. 
* List of transcripts to include in analysis
* Variant info (Name, Position, Strand) as supplied in [docs/ex_variantinfo.txt](docs/ex_variantinfo.txt)

## Installation

Download to your project directory from GitHub, or clone it on the command line using the following command:

```sh
git clone git://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos.git

```

## Information

Script for localizing a given SNV based on its genomic context in the intended CCDS, then generate oligos using the Q5 approach: SNV is in the middle of the F oligo, with at least 10 bases on both sides, and the R oligo is actually on the opposite strand BEFORE the F oligo


			    F oligo -*->
			R oligo <--- 
	 ========================================== vector

The PCR generates a linear version of the vector that can then be ligated.

Simple and efficient.


## Usage

Several files need to be provided in order for the script to work properly.
It requires :
1. An [Exon file](docs/ex_exon_seq.fa) for the transcript of interest. This willl allow to map the variant position from a genomic point of view to a Coding sequence reference. 
2. An [Isoforms file](docs/ex_isoform.txt): List of all possible isoforms of interest for a given gene. 
3. [Variant info](docs/ex_variantinfo.txt): Genomic positioning of each variant, with variant alleles, associated genes, etc. 


Simply calling the script without arguments will use the templates files included in the repository and all default values listed in -help

```sh

python sitedirmutagen.py

```

However, you can specify your Input files of interest by using the various flags

```sh

python sitedirmutagen.py -e path/to/your/exonfile.fa -i path/to/your/isoformfile.txt  -v path/to/your/variantinfofile.txt

```


You can also specify the name of the outputfile you wish to generate for oligos with the -o flag

```sh

python sitedirmutagen.py -o path/to/your/outputfile.txt

```


The script also allows to generate the FASTA sequence for the variants in the variant Info file, for use in sequencing alignment purposes. 
One must first 1. Activate Fasta mode with the -sf flag and then 2. indicate in which directory save those FASTA files with the -vf flag.

```sh

python sitedirmutagen.py -sf -vf path/to/your/fasta/directory/

```


Use the -h flag for further options and help


## In Development



## Known Bugs
IPYTHON has problems with the argparse module, other Python distributions do not. I suggest using Official Python in the meantime. 
Script has not been tested in Python 3.


## Support

Please [open an issue](https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos.git/issues/new) for support.


## Contributing

Please contribute using [Github Flow](https://guides.github.com/introduction/flow/). Create a branch, add commits, and [open a pull request](https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos/compare/).
