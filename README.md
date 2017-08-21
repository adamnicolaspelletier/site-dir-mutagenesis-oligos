# site-dir-mutagenesis-oligos
Generate oligos for site directed mutagenesis with a Q5 Phusion-type home-made approach in an automated fashion 


### Author: Adam-Nicolas Pelletier
### Github: https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos.git

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Information](#information)
- [Usage](#usage)
- [In Development](#in-development)
- [Known Bugs](#known-bugs)
- [Support](#support)
- [Contributing](#contributing)



## Requirements

* Python Libraries :
	1. Pandas (+ Numpy) 0.20.1+
	2. BioPython 1.69+
	3. adamP_BioTools 1.16+

* Exon FASTA sequences for genes or transcripts of interest. 
* List of transcripts to include in analysis
* Variant info (Name, Position, Strand) as supplied in docs/ex_snvinfo.txt



## Installation

Download to your project directory from GitHub, or clone it on the command line using the following command:

```sh
git clone git://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos.git

```

## Information

Script for localizing a given SNV based on its genomic context in the intended CCDS, then generate oligos using the Q5 approach: SNV is in the middle of the F oligo, with at least 10 bases on both sides, and the R oligo is actually on the opposite strand BEFORE the F oligo


			    F oligo -*->
			R oligo <--- 

vector: ====================================================================

The PCR generates a linear version of the vector that can then be ligated.

Simple and efficient.



## Usage

This script was designed with high-throughput in mind: it will generate a high amount of oligos for Site-directed mutagenesis ina short time. 
However, to do so,it requires information: 

1. The fasta sequence of the [exon reference](docs/ex_exon_seq.fa) sequences can be specified using the -e flag or --exonfile. 

```sh

	python sitedirmutagen.py -e path/to/my/exonfastafile.fa -f 
com
```


2. The list of possible [isoforms](docs/ex_ex_isoforms.txt). This allows you to choose which isoform to mutate. 


```sh

	python sitedirmutagen.py -i path/to/my/isoform file

```
3. [Variant information](docs/ex_variantinfo.txt). This will contain position, variant switches, strand, etc. 


```sh

	python sitedirmutagen.py -snv path/to/my/isnvinfofile

```

Additional parameters are also available to specify output file name, oligo length,  and variant fasta file directory using -o , -ol and -sf respectively.
Please include a trailing / after the end of the variant fasta file directory

```sh

	python sitedirmutagen.py -o path/to/my/outputfile -ol 21 -sf path/to/my/variantfasta/output/directory/

```


Use the -h flag for further options and help



## In Development
Nothing in development at this moment


## Known Bugs
No bugs known


## Support

Please [open an issue](https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos.git/issues/new) for support.


## Contributing

Please contribute using [Github Flow](https://guides.github.com/introduction/flow/). Create a branch, add commits, and [open a pull request](https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos/compare/).



BUGS:
If you find a bug or have a suggestion please end your informations to adam.nicolas.pelletier@gmail.com