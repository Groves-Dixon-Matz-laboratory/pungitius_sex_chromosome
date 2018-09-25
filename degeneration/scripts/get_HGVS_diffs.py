#!/usr/bin/env python
##get_HGVS_diffs.py
##written 8/16/18
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
from Bio import AlignIO




##############################
###### DEFINE FUNCTIONS ######
##############################



##################################
############## MAIN ##############
##################################


if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Takes a multiple alingment and two seuqnce IDs from that alignment.
    Outputs the differences between them in Human Genome Variation Society format eg:
    Q5A = Q at position 5 is changed to A

    - symbol is considered a deletion


    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input_file', help = 'The the input file')
    parser.add_argument('-s1', required = True, dest = 'alignment_id1', help = 'The subject sequence in alignment you want to compare. Think of this as reference.')
    parser.add_argument('-s2', required = True, dest = 'alignment_id2', help = 'The query sequence. If s1 has G and s2 has A at site 4, output will be G4A')
    parser.add_argument('-v', required = False, default = False, dest = 'verbose', help = 'Print lots of stuff to standart out while running.')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infileName = args.input_file
    s1 = args.alignment_id1
    s2 = args.alignment_id2
    verbose = args.verbose


    alignment = AlignIO.read(infileName, "fasta")
    for record in alignment:
        if record.id==s1:
            sr1 = record
        if record.id==s2:
            sr2 = record
    if verbose:
        print("\nRecord for s1:")
        print(sr1)
        print("\nRecord for s2:")
        print(sr2)

    str1 = str(sr1.seq)
    str2 = str(sr2.seq)

    if not len(str1) == len(str2):
        exit("Error. Sequence lengths are not equal")

    matchCount = 0
    missCount = 0
    missList = []
    for i in range(len(str1)):
        aa1 = str1[i]
        aa2 = str2[i]
        if aa1==aa2:
            matchCount += 1
            continue
        else:
            missCount += 1
            pos = i+1
            hgvs = "{}{}{}".format(aa1, pos, aa2)
            missList.append(hgvs)
    if verbose:
        print("\nAlignment lengths = {}".format(len(str1)))
        print("Number matching   =  {}".format(matchCount))
        print("Number mismatch   =  {}".format(missCount))
        if (missCount + matchCount) != len(str1):
            exit("Error. Missing some calls :(")
        else:
            print("All sites accounted for.")
        print("Missmatches:")
        for mm in missList:
            print(mm)
            
    outName = "{}_{}_HGVS_diffs_{}.txt".format(s1, s2, infileName.rstrip(".fasta"))
    if len(missList) > 0:
        if verbose:
            print("\nSaving as {}".format(outName))
        with open(outName, 'w') as out:
            for mm in missList:
                out.write(mm + "\n")
    else:
        if verbose:
            print("No variants found. Skipping output.")
    if verbose:
        print("\nDone.\n")













