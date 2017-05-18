  

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### scriptname.py
### Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 

import numpy as np
import itertools
import random
import os
import pandas as pd
from Bio import SeqIO
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline


##Goal: to take snp data and design oligos for site directed mutagenesis on WT vectors. 

########################################################################################################################################################
########################################################## USER INPUT ##################################################################################
snpinfo = "snpinfo.txt"
isoformfile = "isoforms.txt"

letters = {0:"A",1:"C",2:"G",3:"T"}  # dictionary of indexes of each nucleotide for matrices
lettersinv = {"A":0,"C":1,"G":2,"T":3}
genefile = "gene_seq.fasta"
cdnafile = "cdna_seq.fasta"
ccdsfile = "ccds_seq.fasta"
exonfile = "exon_seq.fasta"
primeroutput = "SDM_primer_output.txt"

########################################################################################################################################################
########################################################################################################################################################

#

def fastaconvert(fastalist,removepos):  
    """convert a conventional fasta file into a list of IDs and whole sequences (merges the 50 characters per line)
    takes a listfrom readlines. the removepos argument allows the user to remove items from the list, eithher the first 1 (0) , the last 1 (-1) , or none (NO) """
    a = ""
    x = ""
    z = []
    for i in fastalist:
        if ">" in i:
            a = a.replace("\n","")  
            z.append(a)
            x += i
            x = x.replace("\n","")
            z.append(x)
            a = ""
            x = ""
        else:
            a += i
    if removepos == "NO":
        pass
    else: 
        del z[removepos]
    return z

def randomseq(length, format):   
    """ generates a random n length DNA sequence as a numpy array"""
    matrix = np.zeros( (4, length) )
    index = []
    for i in range (length):    
        index.append([random.randrange(0,4,1), i]) 
    a = np.array(index)  
    matrix[a[:,0], a[:,1]] = 1
    if format == "numpy":
        return matrix
    elif format == "string":
        return matrixmaker(matrix)

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

# promoterlist = []
# for i in os.listdir("promoter_files/seeded/"+str(n)+"bp/"):
#     if i.endswith(".txt"):
#         promoterlist.append("promoter_files/seeded/"+str(n)+"bp/"+i)


fasta = open("tfcds.txt")  #name of the input fasta one liner file. 
fastalist = fasta.readlines()
fastalist.append(">")

fastagenes = open("tfgenes.txt")
fastalistgenes = fastagenes.readlines()
fastalistgenes.append(">")


fastafinal = fastaconvert(fastalist, -1)
fastafinalgenes = fastaconvert(fastalistgenes, -1)

snpinfodf = pd.read_csv(snpinfo, sep = "\t")
exondb = pd.read_csv("ExonDB.txt", sep = "\t")
titles = ["EnsGeneID","EnsTransID","Exonstart","Exonstop","ExonID","Rank_in_transcript","Sequence"]

isoforms = open(isoformfile).read().splitlines()
snpinfodf2 = snpinfodf[snpinfodf["Ensembl Transcript ID"].isin(isoforms)]
snps = snpinfodf2["Variation Name"].tolist()


snpoutputint = snpinfodf2
dfindex = range(len(snpoutputint))
snpoutputint["index"] = dfindex
snpoutput = snpoutputint.set_index("index")
# print snpoutput


# print snpoutputint
# print snpoutput



fwdprimer = []
revprimer = []
sequenced = []



print snps
print 


progress = 0.0
totalprog = float(len(snpoutput))
print len(snpoutput)

# with open("SNVfasta.txt", "w") as snvfasta:
#     snvfasta.write("")
with open("tempoutput.txt", "w") as tempoutput:
    tempoutput.write("")
for i in xrange(len(snps)):

    progress += 1
    
    fastadict = {}
    exontotal= []
    for recorda in SeqIO.parse(exonfile, "fasta"):
            if snpinfodf2.iloc[i]["Ensembl Gene ID"] in recorda.id:
                idsplit = (str(recorda.id).split("|"))
                idsplit.append(str(recorda.seq))
                idsplit[2] = int(idsplit[2])
                idsplit[3] = int(idsplit[3])
                exontotal.append(idsplit)
                          
                exondb = pd.DataFrame(exontotal, columns=titles)
    # print exondb
    
               
    try:
        start = int(snpinfodf2.iloc[i]["Gene_start"])
        snppos =  int(snpinfodf2.iloc[i]["Position on Chromosome (bp)"])
        exonsnpdf = exondb[(exondb["Exonstart"]<snppos) & (exondb["Exonstop"]>snppos)]
        
        # print snppos
        # print exonsnpdf
        genename = snpinfodf2.iloc[i]["GENENAME"]  
        insertfile = "INSERTS/%s.txt" %genename

        


        original = ""
        mutantseq = ""

        for j in xrange(len(exondb)):
            exseq = (exondb[exondb["Rank_in_transcript"] == str(j+1)])
            original += str(exseq.iloc[0]["Sequence"])


        if snpinfodf2.iloc[i]["Gene_strand"] == -1:
            oriseq = list(reversecomp((exonsnpdf.iloc[0]["Sequence"])))
            mutlist = oriseq
            # print oriseq[(snppos-exonsnpdf.iloc[0]["Exonstart"])]
            # print mutlist[(snppos-exonsnpdf.iloc[0]["Exonstart"])]
            mutlist[(snppos-exonsnpdf.iloc[0]["Exonstart"])] = snpinfodf2.iloc[i]["Variant_allele"]
            mutantseq = reversecomp("".join(mutlist))

            mutant = original.replace(exonsnpdf.iloc[0]["Sequence"], mutantseq)


            
            # print exondb
            # print reversecomp((exonsnpdf.iloc[0]["Sequence"]))[(snppos-exonsnpdf.iloc[0]["Exonstart"])] , snpinfodf2.iloc[i]["Classis_allele"]
            # print snpinfodf2.iloc[i]["Classis_allele"], snpinfodf2.iloc[i]["Variant_allele"]
            # print "-"

        else:
            mutlist= list((exonsnpdf.iloc[0]["Sequence"]))
            mutlist[(snppos-exonsnpdf.iloc[0]["Exonstart"])] = snpinfodf2.iloc[i]["Variant_allele"]
            mutantseq = "".join(mutlist)
 
            mutant = original.replace((exonsnpdf.iloc[0]["Sequence"]),mutantseq)
            # print snpinfodf2.iloc[i]["Classis_allele"], snpinfodf2.iloc[i]["Variant_allele"]
            # print "+"
       

        for z in xrange(len(original)):
            if list(original)[z] != list(mutant)[z]:
                # print z
                # print list(original)[z], list(mutant)[z]
                fwdprimer.append("".join(mutant[z-10:z+11]))
                revprimer.append(reversecomp("".join(mutant[z-31:z-10])))
               
            
                # print revprimerls
                # print "".join(mutant[z-31:z+11])

                
            else:
                pass

        try:
            insertfile = "INSERTS/%s.txt" % genename
            for recordc in SeqIO.parse(insertfile, "fasta"):
                a = 1
            sequenced.append("YES")

        except IOError:
            sequenced.append("NO")
        
        # with open("tempoutput.txt", "a") as tempoutput:
        #     tempoutput.write(">1. Original  "+ str(recorda.id)+ "\n" + original +"\n")
        with open("tempoutput.txt", "a") as tempoutput:
            tempoutput.write(">"+ str(snpinfodf2.iloc[i]["Variation Name"])+ "\n" + mutant +"\n")
        print str(progress) , snpinfodf2.iloc[i]["Variation Name"], "DONE", str(progress/totalprog*100) , "%"


        

    except IndexError:
        print str(progress) , snpinfodf2.iloc[i]["Variation Name"], "NOT FOUND INSIDE CODING EXON OF CHOSEN TRANSCRIPT" , str(progress/totalprog*100) , "%"
        fwdprimer.append("---")
        revprimer.append("---")
        sequenced.append("NO")
        pass
print snpinfodf2

snpoutput["FWD_Primer"] = fwdprimer
snpoutput["REV_Primer"] = revprimer
snpoutput["SEQUENCED"] = sequenced
snpoutput.to_csv(primeroutput, sep='\t')
print snpoutput




   