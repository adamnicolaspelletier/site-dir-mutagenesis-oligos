site-dir-mutagenesis-oligos

Generate oligos for site directed mutagenesis with a Q5 Phusion-type home-made approach in an automated fashion

Script for localizing a given SNV based on its genomic context in the intended CCDS, then generate oligos using the Q5 approach: SNV is in the middle of the F oligo, with at least 10 bases on both sides, and the R oligo is actually on the opposite strand BEFORE the F oligo


			    F oligo -*->
			R oligo <--- 

vector: ====================================================================

The PCR generates a linear version of the vector that can then be ligated.

Simple and efficient.


Requirements:
  - Python Libraries :

	Pandas (+ Numpy)
	BioPython

 - Exon FASTA sequences for genes or transcripts of interest. 
 - List of transcripts to include in analysis
 - Variant info (Name, Position, Strand) as supplied in docs/ex_snvinfo.txt



Use the -h flag for further options and help




BUGS:
If you find a bug or have a suggestion please end your informations to adam.nicolas.pelletier@gmail.com