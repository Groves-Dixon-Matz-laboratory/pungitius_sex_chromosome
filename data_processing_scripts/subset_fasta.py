#!/usr/bin/env python
#subset_fasta.py
##written 7/31/15 by Groves Dixon

#Import modules
import time
from sys import argv
from sys import exit
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


usage = "subset_fasta.py fileName.fasta subset_of_sequences.txt > subset_output.fasta"



try:
    if argv[1] == '-h':
        exit("\nUsage:\n" + usage)
except IndexError:
    exit("\nUsage:\n" + usage)

fasta = argv[1]   #the fasta file your subsetting
subset = argv[2]  #single column text table of the sequence names you want to keep



subsetList = []
with open(subset, 'r') as infile:
    for line in infile:
        subsetList.append(line.strip("\n"))
fasSeqs = SeqIO.parse(open(fasta), 'fasta')
for i in fasSeqs:
    if i.id in subsetList:
        print ">" + i.id
        print i.seq
    else: continue

       

