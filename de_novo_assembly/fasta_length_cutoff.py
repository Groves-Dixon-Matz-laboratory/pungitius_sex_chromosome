#!/usr/bin/env python
##fasta_length_cutoff.py


Description = '''
Description:
Returns fasta sequences only above a given threshold length
'''

AdditionalProgramInfo = '''
Additional Program Information:

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo)
parser.add_argument('-fa', required = True, dest = 'fa', help = 'The fasta file')
parser.add_argument('-cut', required = True, dest = 'cutoff', help = 'Length cutoff')
parser.add_argument('-o', required = True, dest = 'outfileName', help = 'Name for output')
args = parser.parse_args()

#Assign Arguments
fastaFile = args.fa
cut = int(args.cutoff)
outfileName = args.outfileName

#iterate through the seqs
fasSeqs = SeqIO.parse(open(fastaFile), 'fasta')
shortRemoved = 0
totalSeqs = 0
kept = 0
print("\nReading through fasta file {}...".format(fastaFile))
with open(outfileName, 'w') as out:
    for seq in fasSeqs:
        totalSeqs += 1
        length = len(seq.seq)
        if length >= cut:
            kept += 1
            out.write(">" + seq.id + "\n" + str(seq.seq) + "\n")
        else:
            shortRemoved += 1
            continue

print("\n\n{} total sequences found in file.".format(totalSeqs))
print("\n{} sequences were greater than or equal to {} in length and kept.".format(kept, cut))
print("\n{} sequences were less than {} and were removed".format(shortRemoved, cut))
print("\nRetained sequences saved as {}".format(outfileName))  
    
