  

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### sitedirmutagen.py
### Copyright 2017 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 2017-05-24

import numpy as np
import os
import pandas as pd
from Bio import SeqIO
import argparse
import time 

pd.options.mode.chained_assignment = None  # default='warn

##Goal: to take variant data and design oligos for site directed mutagenesis on WT vectors. 




parser = argparse.ArgumentParser(description="""Generate oligos for site directed mutagenesis with a Q5 Phusion-type home-made approach in an automated fashion

Script for localizing a given variant based on its genomic context in the intended CCDS, then generate oligos using the Q5 approach: 
variant is in the middle of the F oligo, with at least 10 bases on both sides, and the R oligo is actually on the opposite strand BEFORE the F oligo

Supply FASTAS files for Genes and CDS (Coding Sequence), the variant info as shown in the variantinfo.txt file, """ )

## !!Add metadata in header for filenames, date, command used (" python sitedirmutagene.py -fg "fastagene.fa" -fc fastacds.fa etc")

cwd = os.path.dirname(os.path.realpath(__file__))

parser.add_argument("-e","--exonfile",
                    help="File Containing reference Exon FASTA sequences. Can be obtained from BioMart, or other sources. Defaults to 'docs/ex_exon_seq.fa'", default= "docs/ex_exon_seq.fa")
parser.add_argument("-i", "--isoform", default=str(cwd)+"docs/ex_isoform.txt",
                    help="List of possible isoforms file. Defaults to 'docs/ex_isoform.txt'")
parser.add_argument("-v", "--variantinfo", default=str(cwd)+"docs/ex_variantinfo.txt",
                    help="File containing variant info (Name, Position, Variants), as shown in docs/variantinfo.txt. Can be obtained from BioMart, or other sources. Defaults to 'docs/ex_variantinfo.txt'")
parser.add_argument("-o", "--output", default=str(cwd)+"docs/site_dir_mutagen_output.txt",
                    help="Outputfile containing Oligos for each specified variant within the context of the CDS. Defaults to 'docs/site_dir_mutagen_output.txt'")
parser.add_argument("-ol", "--oligolen", default=21,
                    help="Length of Forward and Reverse Oligos for SDM. Defaults to 21")
parser.add_argument("-vf", "--variantfasta", default=str(cwd)+"docs/VARIANT_FASTA/",
                    help="Directory for Variant CDS FASTA files per gene. Requires -sf flag. Will create directory if it does not exist. Defaults to docs/VARIANT_FASTA/")
parser.add_argument("-sf", "--savefasta", action="store_true",
                    help="Activates save fasta mode, to allow the script to output the FASTA sequence of the difference variants to the -vf flag directory")
args = parser.parse_args()


exonfile = args.exonfile
isoform = args.isoform
variantinfo = args.variantinfo
output = args.output
oligolen = args.oligolen
variantfasta = args.variantfasta
savefasta = args.savefasta


print "\nsitedirectedmutagen.py script for automating mutagenesis oligo generating in high-throughput experiments\n"
print " ** use the -h flag for detailed help on how to use this script **\n"


print "\n\nUsing %s for Reference Exon FASTA ..." % exonfile
print "Using %s for Possible Isoforms File ..." % isoform
print "Using %s for Additional Variant Information ..." % variantinfo
print "Using %s for Output file ..." % output
print "Using %s for Oligo Length ..." % oligolen
if savefasta == True:
    print "Using %s directory for Variant FASTA files  ... " % variantfasta
else:
    print "Not Saving Variant CDS FASTA sequences..."

print "\n\n"



def reversecomp(rprimsequence): ## make a complement version of the sequence, and reverse it so it has the proper orientation
    a = ""
    tempzrev = rprimsequence
    tempzrev = tempzrev.replace("T","X")
    tempzrev = tempzrev.replace("A","T")
    tempzrev = tempzrev.replace("X","A")
    tempzrev = tempzrev.replace("C","Y")
    tempzrev = tempzrev.replace("G","C")
    tempzrev = tempzrev.replace("Y","G")
    templist = list(tempzrev)
    templist.reverse()
    for i in templist:
        a += i
    return a



variantinfodf = pd.read_csv(variantinfo, sep = "\t")    ## Read variantinfo in pandas dataframe
isoforms = open(isoform).read().splitlines()  ## Make list of isoforms
variantinfodf2 = variantinfodf[variantinfodf["Ensembl Transcript ID"].isin(isoforms)]  #Filter variant info based on selected isoforms


variantinfodf2["index"] = range(len(variantinfodf2))
variantinfodf2 = variantinfodf2.set_index("index")  #reindex after filtering, to facilitate .iloc indexing


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
    os.makedirs(variantfasta)
except OSError:
    if not os.path.isdir(variantfasta):
        raise 

if savefasta == True:
    genelist = list(variantinfodf2["GENENAME"].unique())
    for i in genelist:  
        insertfile = "%s%s.txt" % (variantfasta,i)
        with open(insertfile, "w") as variantfastaoutput:
            variantfastaoutput.write("")
        


## Progress Count for oligo design
progress = 1.0
totalprog = float(len(variantinfodf2))

### Lists of Forward, Reverse and 
fwdprimer = []
revprimer = []


### 
for i in xrange(len(variantinfodf2["Variation Name"])):  #Iterate though all variants in variant DataFrame. 
    
    try:
        # start = int(variantinfodf2.iloc[i]["Gene_start"])
        geneid = variantinfodf2.iloc[i]["Ensembl Gene ID"]
        variantpos =  int(variantinfodf2.iloc[i]["Position on Chromosome (bp)"])

        exondbfilt = exondb[exondb["EnsGeneID"]==geneid] #filter Exon dataframe to the gene associated with the variant in the current iteration
        exonvariantdf = exondbfilt[(exondbfilt["Exonstart"]<variantpos) & (exondbfilt["Exonstop"]>variantpos)] # Filter the Exon DF to the Exon that fits the genome coordinates
        

        original = ""
        mutantseq = ""

        for j in xrange(len(exondbfilt)):  #reconstitute CDS from sequence of exons in order
            exseq = (exondbfilt[exondbfilt["Rank_in_transcript"] == str(j+1)])
            original += str(exseq.iloc[0]["Sequence"])


        if variantinfodf2.iloc[i]["Gene_strand"] == -1: ## Variants are always on the + strand: if gene is on the - strand, reverse complement genesequence to fit variant properly. 
            mutlist = list(reversecomp((exonvariantdf.iloc[0]["Sequence"])))
            mutlist[(variantpos-exonvariantdf.iloc[0]["Exonstart"])] = variantinfodf2.iloc[i]["Variant_allele"] # replace 
            mutantseq = reversecomp("".join(mutlist)) #mutated exon
            mutant = original.replace(exonvariantdf.iloc[0]["Sequence"], mutantseq) # replace the mutated exon in the original sequence
        else:
            mutlist= list((exonvariantdf.iloc[0]["Sequence"]))
            mutlist[(variantpos-exonvariantdf.iloc[0]["Exonstart"])] = variantinfodf2.iloc[i]["Variant_allele"]
            mutantseq = "".join(mutlist) #mutated exon
            mutant = original.replace((exonvariantdf.iloc[0]["Sequence"]),mutantseq) # replace the mutated exon in the original sequence


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
        


        genename = variantinfodf2.iloc[i]["GENENAME"]
        ### if User chose to output the Variant FASTA sequences to specified directory only
        if savefasta == True:
            try:
                currdate = time.strftime("%d/%m/%Y")
                variantfastafile = "%s%s.txt" % (variantfasta,genename)
                if os.stat(variantfastafile).st_size == 0: # if specified file is empty add Reference Consensus sequence
                    out = ">"+ str(variantinfodf2.iloc[i]["Ensembl Gene ID"])+ "|"+ str(variantinfodf2.iloc[i]["Ensembl Transcript ID"])+ "|" + \
                    str(variantinfodf2.iloc[i]["GENENAME"])+ "| REFERENCE | Date: " + str(currdate) + "|"+ "sitedirmutagen.py" + "|"+ variantinfo + "|"+ exonfile + "|" + \
                    isoform + "\n" + original +"\n"
                    with open(variantfastafile, "a") as tempoutput:
                        tempoutput.write(out)

                with open(variantfastafile, "a") as tempoutput: # append variant fasta sequences
                    out = ">"+ str(variantinfodf2.iloc[i]["Variation Name"])+ "|"+ str(variantinfodf2.iloc[i]["Ensembl Gene ID"]) + "|"+ \
                    str(variantinfodf2.iloc[i]["Ensembl Transcript ID"]) + "|"+ str(variantinfodf2.iloc[i]["GENENAME"]) + "\n" + mutant +"\n"
                    tempoutput.write(out)

            except IOError:
                sequenced.append("NO")

        print str(progress) , variantinfodf2.iloc[i]["Variation Name"], "DONE", str(progress/totalprog*100) , "%" ## Progress report on total variant list


    except IndexError:
        print str(progress) , variantinfodf2.iloc[i]["Variation Name"], "NOT FOUND INSIDE CODING EXON OF CHOSEN TRANSCRIPT" , str(progress/totalprog*100) , "%"
        fwdprimer.append("---")
        revprimer.append("---")
    progress += 1
        


variantinfodf2["FWD_Primer"] = fwdprimer
variantinfodf2["REV_Primer"] = revprimer
variantinfodf2.to_csv(output, sep='\t')  # add oligos to dataframe and output the dataframe to specified file 
