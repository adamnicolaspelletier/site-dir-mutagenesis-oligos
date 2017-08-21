  

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### sitedirmutagen.py
### Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 2017-05-24

import numpy as np
import os
import pandas as pd
from Bio import SeqIO
import argparse
import time 
from adamP_BioTools.dna_tools import *
from adamP_BioTools.os_handling import path_validity

pd.options.mode.chained_assignment = None  # default='warn

##Goal: to take snv data and design oligos for site directed mutagenesis on WT vectors. 



########################################################################################################################################################
########################################################## USER INPUT ##################################################################################

parser = argparse.ArgumentParser(description="""Generate oligos for site directed mutagenesis with a Q5 Phusion-type home-made approach in an automated fashion

Script for localizing a given SNV based on its genomic context in the intended CCDS, then generate oligos using the Q5 approach: 
SNV is in the middle of the F oligo, with at least 10 bases on both sides, and the R oligo is actually on the opposite strand BEFORE the F oligo

Supply FASTAS files for Genes and CDS (Coding Sequence), the variant info as shown in the snvinfo.txt file, """ )

## !!Add metadata in header for filenames, date, command used (" python sitedirmutagene.py -fg "fastagene.fa" -fc fastacds.fa etc")

parser.add_argument("-e","--exonfile",
                    help="File Containing reference Exon FASTA sequences. Can be obtained from BioMart, or other sources. Defaults to 'docs/ex_exon_seq.fa'", default= "docs/ex_exon_seq.fa")
parser.add_argument("-i", "--isoform", default="docs/ex_isoform.txt",
                    help="List of possible isoforms file. Defaults to 'docs/ex_isoform.txt'")
parser.add_argument("-snv", "--snvinfo", default="docs/ex_snvinfo.txt",
                    help="File containing SNV info (Name, Position, Variants), as shown in docs/snvinfo.txt. Can be obtained from BioMart, or other sources. Defaults to 'docs/ex_snvinfo.txt'")
parser.add_argument("-o", "--output", default="docs/site_dir_mutagen_output.txt",
                    help="Outputfile containing Oligos for each specified SNV within the context of the CDS. Defaults to 'docs/site_dir_mutagen_output.txt'")
parser.add_argument("-ol", "--oligolen", default=21,
                    help="Length of Forward and Reverse Oligos for SDM. Defaults to 21")
parser.add_argument("-sf", "--snvfasta", default="docs/VARIANT_FASTA/",
                    help="Directory for Variant CDS FASTA files per gene. Requires 'y' when user is prompted. Will create dircetory if it does not exist. Defaults to docs/VARIANT_FASTA/")

args = parser.parse_args()


exonfile = args.exonfile
isoform = args.isoform
snvinfo = args.snvinfo
output = args.output
oligolen = args.oligolen
snvfasta = args.snvfasta

print "\nsitedirectedmutagen.py script for automating mutagenesis oligo generating in high-throughput experiments\n"
print " ** use the -h flag for detailed help on how to use this script **\n"


# Prompt user for exporting Variant CDS FASTA sequences to the VARIANT_FASTA directory
varfasta = raw_input("Save Variant CDS FASTA sequences?  (y/n):  ")


print "\n\nUsing %s for Reference Exon FASTA ..." % exonfile
print "Using %s for Possible Isoforms File ..." % isoform
print "Using %s for Additional SNV Information ..." % snvinfo
print "Using %s for Output file ..." % output
print "Using %s for Oligo Length ..." % oligolen
if varfasta.upper() == "Y":
    print "Using %s directory for Variant FASTA files  ... " % snvfasta
else:
    print "Not Saving Variant CDS FASTA sequences..."

print "\n\n"

########################################################################################################################################################
########################################################################################################################################################


#file verification step. Although not necessary in most cases, it is frustrating to get an OSError for trying to export the output to a file halfway through a long loop.
if path_validity(exonfile) == False:
    raise OSError("Invalid file or directory: %s" % exonfile) 
elif path_validity(isoform) == False:
    raise OSError("Invalid file or directory: %s" % isoform) 
elif path_validity(snvinfo) == False:
    raise OSError("Invalid file or directory: %s" % snvinfo) 
elif path_validity(output) == False:
    raise OSError("Invalid file or directory: %s" % output) 





snvinfodf = pd.read_csv(snvinfo, sep = "\t")    ## Read SNVinfo in pandas dataframe
isoforms = open(isoform).read().splitlines()  ## Make list of isoforms
snvinfodf2 = snvinfodf[snvinfodf["Ensembl Transcript ID"].isin(isoforms)]  #Filter SNV info based on selected isoforms


snvinfodf2["index"] = range(len(snvinfodf2))
snvinfodf2 = snvinfodf2.set_index("index")  #reindex after filtering, to facilitate .iloc indexing


#Create  temporary file to output preliminary data. 
with open("tempoutput.txt", "w") as tempoutput:
    tempoutput.write("")



### Make Exon DataFrame from exon file .
titles = ["EnsGeneID","EnsTransID","Exonstart","Exonstop","ExonID","Rank_in_transcript","Sequence"]
exontotal= []
for record in SeqIO.parse(exonfile, "fasta"):
    idsplit = (str(record.id).split("|")) # Split GeneID, Transcript ID, Exon Genomic start pos and Exon Genomic Stop Pos. 
    idsplit.append(str(record.seq)) # Add FASTA sequence to this list
    idsplit[2] = int(idsplit[2])
    idsplit[3] = int(idsplit[3])
    exontotal.append(idsplit) # List of lists with all exon info
                  
exondb = pd.DataFrame(exontotal, columns=titles)

try:
    os.makedirs(snvfasta)
except OSError:
    if not os.path.isdir(snvfasta):
        raise 

if varfasta.upper() == "Y":
    genelist = list(snvinfodf2["GENENAME"].unique())
    for i in genelist:  
        insertfile = "%s%s.txt" % (snvfasta,i)
        with open(insertfile, "w") as snvfastaoutput:
            snvfastaoutput.write("")
        


## Progress Count for oligo design
progress = 1.0
totalprog = float(len(snvinfodf2))

### Lists of Forward, Reverse and 
fwdprimer = []
revprimer = []


### 
for i in xrange(len(snvinfodf2["Variation Name"])):  #Iterate though all SNVs in SNV DataFrame. 
    
    try:
        # start = int(snvinfodf2.iloc[i]["Gene_start"])
        geneid = snvinfodf2.iloc[i]["Ensembl Gene ID"]
        snvpos =  int(snvinfodf2.iloc[i]["Position on Chromosome (bp)"])

        exondbfilt = exondb[exondb["EnsGeneID"]==geneid] #filter Exon dataframe to the gene associated with the SNV in the current iteration
        exonsnvdf = exondbfilt[(exondbfilt["Exonstart"]<snvpos) & (exondbfilt["Exonstop"]>snvpos)] # Filter the Exon DF to the Exon that fits the genome coordinates
        

        original = ""
        mutantseq = ""

        for j in xrange(len(exondbfilt)):  #reconstitute CDS from sequence of exons in order
            exseq = (exondbfilt[exondbfilt["Rank_in_transcript"] == str(j+1)])
            original += str(exseq.iloc[0]["Sequence"])


        if snvinfodf2.iloc[i]["Gene_strand"] == -1: ## Variants are always on the + strand: if gene is on the - strand, reverse complement genesequence to fit SNV properly. 
            mutlist = list(reversecomp((exonsnvdf.iloc[0]["Sequence"])))
            mutlist[(snvpos-exonsnvdf.iloc[0]["Exonstart"])] = snvinfodf2.iloc[i]["Variant_allele"] # replace 
            mutantseq = reversecomp("".join(mutlist)) #mutated exon
            mutant = original.replace(exonsnvdf.iloc[0]["Sequence"], mutantseq) # replace the mutated exon in the original sequence
        else:
            mutlist= list((exonsnvdf.iloc[0]["Sequence"]))
            mutlist[(snvpos-exonsnvdf.iloc[0]["Exonstart"])] = snvinfodf2.iloc[i]["Variant_allele"]
            mutantseq = "".join(mutlist) #mutated exon
            mutant = original.replace((exonsnvdf.iloc[0]["Sequence"]),mutantseq) # replace the mutated exon in the original sequence


        if oligolen % 2 == 0:
            flankf = (oligolen/2)  #if length of oligo is even, add same number of bases on both sides of the variant. if not, add 1 base on the 3' side
        else:
            flankf = (oligolen/2) +1
        
        flankr = (oligolen/2)


        for z in xrange(len(original)): #Iterate through the sequence until the variant base is found.
            if list(original)[z] != list(mutant)[z]:                      
                fwdprimer.append("".join(mutant[z-flankr:z+flankf])) ## F primer has variant base in middle, and flanked by oligolength /2
                revprimer.append(reversecomp("".join(mutant[z-(flankr+oligolen):z-flankr]))) ### R primer begins in 3' where the 5' portion of F primer begins, no overlap. 
            else:
                pass
        


        genename = snvinfodf2.iloc[i]["GENENAME"]
        ### if User chose to output the Variant FASTA sequences to specified directory only
        if varfasta.upper() == "Y":
            try:
                currdate = time.strftime("%d/%m/%Y")
                snvfastafile = "%s%s.txt" % (snvfasta,genename)
                if os.stat(snvfastafile).st_size == 0: # if specified file is empty add Reference Consensus sequence
                    out = ">"+ str(snvinfodf2.iloc[i]["Ensembl Gene ID"])+ "|"+ str(snvinfodf2.iloc[i]["Ensembl Transcript ID"])+ "|" + \
                    str(snvinfodf2.iloc[i]["GENENAME"])+ "|"+ "Date: " + str(currdate) + "|"+ "sitedirmutagen.py" + "|"+ snvinfo + "|"+ exonfile + "|" + \
                    isoform + "\n" + original +"\n"
                    with open(snvfastafile, "a") as tempoutput:
                        tempoutput.write(out)

                with open(snvfastafile, "a") as tempoutput: # append variant fasta sequences
                    tempoutput.write(">"+ str(snvinfodf2.iloc[i]["Variation Name"])+ "\n" + mutant +"\n")

            except IOError:
                sequenced.append("NO")

        print str(progress) , snvinfodf2.iloc[i]["Variation Name"], "DONE", str(progress/totalprog*100) , "%" ## Progress report on total SNV list


    except IndexError:
        print str(progress) , snvinfodf2.iloc[i]["Variation Name"], "NOT FOUND INSIDE CODING EXON OF CHOSEN TRANSCRIPT" , str(progress/totalprog*100) , "%"
        fwdprimer.append("---")
        revprimer.append("---")
    progress += 1
        


snvinfodf2["FWD_Primer"] = fwdprimer
snvinfodf2["REV_Primer"] = revprimer
snvinfodf2.to_csv(output, sep='\t')  # add oligos to dataframe and output the dataframe to specified file 
