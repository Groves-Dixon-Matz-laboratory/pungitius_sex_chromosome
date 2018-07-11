#!/usr/bin/env python
##fasta_sequence_characters.py
##written 11/4/15 by Groves Dixon
ProgramName = 'fasta_sequence_characters.py'
LastUpdated = '11/4/15'
By = 'Groves Dixon'
VersionNumber = '1.0'

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:

Returns information about the sequences in a fasta.
Currently just returns length.


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
args = parser.parse_args()

#Assign Arguments
Fa = args.fa

#print the header
# print "contig\tlength"  

#iterate through the seqs
fasSeqs = SeqIO.parse(open(Fa), 'fasta')
for seq in fasSeqs:
    print "%s\t%i" % (seq.id, len(seq.seq))
    cgCount = str(seq.seq).count('CG')
    # print "CpG count = {}".format(cgCount)
    


#return time to run
# print('\nTime took to run: {}'.format(Time))