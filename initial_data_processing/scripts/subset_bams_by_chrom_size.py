#!/usr/bin/env python

#Get_files
'''
'''

from sys import argv

debug = False

#set up usage information
usage = "fasta_sequence_characters.py -fa $stickleGenome > chromLengths.txt\nsubset_bams_by_chrom_size.py chromLengths.txt 4 *sorted.bam > subsetBams"
try:
    argv[1]
except IndexError:
    print('\nUsage:')
    exit(usage)




infileName = argv[1]
nSubs = int(argv[2])
bams = argv[3:]

if debug: #print to see if list is correct
    print FileList


with open(infileName, 'r') as infile:
	for line in infile:
		line=line.strip("\n").split()
		chrom = line[0]
		length = int(line[1])
		step=length/(nSubs)
		stepList = range(0, nSubs)
		sections = [x*step for x in stepList]
		sections.append(length)
		for b in bams:
			for i in range(len(sections)-1):
				outFile = b.split(".bam")[0] + "_{}_sub{}.bam".format(chrom, i+1) 
				print("samtools view -b {} -o {} {}:{}-{}".format(b, outFile, chrom, sections[i], (sections[i+1]-1)))



