# site-dir-mutagenesis-oligos

### Author: Adam-Nicolas Pelletier
### Github: https://github.com/adamnicolaspelletier/site-dir-mutagenesis-oligos

## Generate oligos for site directed mutagenesis with a Q5 Phusion-type home-made approach in an automated fashion

Script for localizing a given SNV based on its genomic context in the intended CCDS, then generate oligos using the Q5 approach: SNV is in the middle of the F oligo, with at least 10 bases on both sides, and the R oligo is actually on the opposite strand BEFORE the F oligo


			    F oligo -*->
			R oligo <--- 
	 ========================================== vector

The PCR generates a linear version of the vector that can then be ligated.

Simple and efficient.


## Requirements:
  1. Python Libraries :

	* Pandas (+ Numpy)
	* BioPython

 2. Exon FASTA sequences for genes or transcripts of interest. 
 3. List of transcripts to include in analysis
 4. Variant info (Name, Position, Strand) as supplied in docs/ex_snvinfo.txt



Use the -h flag for further options and help



## KNOWN BUGS:
 IPYTHON has problems with the argparse module, other Python distributions do not. I suggest using Official Python in the meantime. 
 Script has not been tested in Python 3.

 If you find a bug or have a suggestion please end your informations to adam.nicolas.pelletier@gmail.com



