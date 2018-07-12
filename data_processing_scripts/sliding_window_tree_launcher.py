#!/usr/bin/env python
##window_vcfkit_output.py
#Groves Dixon
##last revised 1-16-18


#import modules
import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
import math



#main
Start_time = time.time() ##keeps track of how long the script takes to run
##Set Up Argument Parsing
Description = '''
Description:
This script launches sliding window tree-building with RAxML from VCF files by:
1 picking out the windows based on indicated sizes (given as number of SNPs)
2. making commands to pipe VCFs subset to the windows into VCF-kit
3. reading the vcf the fasta came from to output windows' coordinates in long format
4. writing commands to run RAxML on the fasta slices for each window
'''
#designate args
parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = True, dest = 'input', nargs="+", help = 'Glob to VCF files to build windows from')
parser.add_argument('-wSNP', required = False, dest = 'window_size_SNPs', default=False, help = 'Desired window length in SNPs')
parser.add_argument('-wBP', required = False, dest = 'window_size_BPs', default=False, help = 'Desired window length in base-pairs')
parser.add_argument('-s', required = True, dest = 'slide_length', help = 'Slide length. Amount window slides each time either in SNPs or Basepairs depdending on argument for window size')
parser.add_argument('-qTrees', required = False, dest = 'quick_trees', default=False, help = 'Boolean for making quick UPGMA trees with VCF-kit. Default is False')
parser.add_argument('-rTrees', required = False, dest = 'rax_trees', default=False, help = 'Boolean for making raxml trees from the fastas as they are made')
#assign
args = parser.parse_args()
vcfList = args.input
windowSizeSNPs=args.window_size_SNPs
windowSizeBPs=args.window_size_BPs
slideLength = int(args.slide_length)
quickTrees=args.quick_trees
raxTrees = args.rax_trees

#set up window units based on user arguments
if windowSizeSNPs:
    if windowSizeBPs:
        exit("\nErrror. Please choose to give window size either in SNPs or Basepairs...\n")
    else:
        windowSize=int(windowSizeSNPs)
        windowUnit = 'snp'
        print("\nUsing window size set at {} SNPs".format(windowSize))
        print("Slide length = {} SNPs".format(slideLength))
else:
    windowSize=int(windowSizeBPs)
    windowUnit = 'bp'
    print("\nUsing window size set at {} Basepairs".format(windowSize))
    print("Slide length = {} Basepairs".format(slideLength))


#make directory to put the window fasta files in
directory="windowFastas_w{}_s{}_{}".format(windowSize, slideLength, windowUnit)
if not os.path.exists(directory):
    os.makedirs(directory)

#make a directory to put the quick upgma trees in 
if quickTrees:
    print("\nQuick trees argument set to True value, will build trees with vcf-kit")
    treeDir = "quickTrees_w{}_s{}_{}".format(windowSize, slideLength, windowUnit)
    if not os.path.exists(treeDir):
        os.makedirs(treeDir)
else:
    print("\nQuickTrees set to False, will only output fasta files for windows")


if raxTrees:
    print("\nWill build trees using RAXML from fastas as they are created.")

#run functions
#first gather positions of the variants from the VCF
with open('slidingWindowTreeCommands_w{}_s{}.txt'.format(windowSize, slideLength), 'w') as out:
    for vcf in vcfList:
        chrom=vcf.split("_")[0]
        vcfPrefix = vcf.rstrip(".vcf")
        posList=[]
        windowList = []
        print("Reading VCF file {}...".format(vcf))
        with open(vcf, 'r') as infile:
            for line in infile:
                #skip header lines
                if line[0] == "#":
                    continue
                #count up variants
                line=line.split()
                pos=int(line[1])
                posList.append(pos)
            print("\nFound {} variants in the VCF file {}".format(len(posList), vcf))


            #build windows by number of snps
            #this is based simply the number of variants in the VCF
            if windowUnit=='snp':

                #set up windows for SNPs
                tot = len(posList) #total positions in VCF
                lastIndex=tot-1    #index from zero
                runnerUp = lastIndex - windowSize
                lefts=range(0, runnerUp, slideLength)
                rights=[x+(windowSize-1) for x in lefts]
                remainder = lastIndex - rights[-1]
                rights[-1] = lastIndex
                print("From position {} to position {}".format(posList[0], posList[-1]))
                print("Slicing VCF in {} windows each {} SNPs in length with slide size = {}, overlap = {}".format(len(rights), windowSize, slideLength, (windowSize-slideLength)))
                print("{} remaining SNPs will be added into final window".format(remainder))
                print("\nWINDOWS:")
                for i in range(5):
                    print("{}\t{}\tsize={}".format(lefts[i], rights[i], (rights[i]-lefts[i]+1)))
                print("/.../")
                for i in range((len(lefts)-5), len(lefts)):
                    print("{}\t{}\tsize={}".format(lefts[i], rights[i], (rights[i]-lefts[i]+1)))
                print("*Note, these window coordinates refer to the list of variant positions indexed by zero,")
                print("but VCFtools --toBp and --fromBp are inclusive, so both the left and right indexed position are included in each window")

                # print("\nFirst window indices:")
                # print("{}\t{}".format(lefts[0], rights[0]))
                # print("First window positions:")
                # print("{}\t{}".format(posList[lefts[0]], posList[rights[0]]))
                # print("\nLast window indices:")
                # print("{}\t{}".format(lefts[-1], rights[-1]))
                # print("Last window positions:")
                # print("{}\t{}".format(posList[lefts[-1]], posList[rights[-1]]))
                for i in range(len(lefts)):
                    windowNumber='w{}'.format(i)
                    l=lefts[i]
                    r=(rights[i])
                    fromBp = posList[l]
                    toBp = posList[r]
                    if raxTrees:
                        outFasta = "{}/{}_{}.fasta".format(directory, vcfPrefix, windowNumber)
                        treeCommand = "/work/02260/grovesd/stampede2/standard-RAxML/raxmlHPC-PTHREADS -s {} -n {} -m GTRCAT -f d -p 12345 -T 1".format(outFasta, outFasta.rstrip(".fasta").split("/")[1])
                        sliceCommand = "vcftools --vcf {} --chr {} --from-bp {} --to-bp {} --recode -c | vk phylo fasta - > {} && {}\n".format(vcf, chrom, fromBp, toBp, outFasta, treeCommand)
                    else:
                        sliceCommand = "vcftools --vcf {} --chr {} --from-bp {} --to-bp {} --recode -c | vk phylo fasta - > {}/{}_{}.fasta\n".format(vcf, chrom, fromBp, toBp, directory, vcfPrefix, windowNumber)
                    out.write(sliceCommand)
                    if quickTrees:
                        treeCommand = "vcftools --vcf {} --chr {} --from-bp {} --to-bp {} --recode -c | vk phylo tree upgma - > {}/{}_{}.newick\n".format(vcf, chrom, fromBp, toBp, treeDir, vcfPrefix, windowNumber)
                        out.write(treeCommand)
                    windowList.append([windowNumber, str(fromBp), str(toBp)])

            #build windows by physical distance
            else:
                print("\nDividing chromosomes into windows based on physical distance")
                print("window size = {}bp".format(windowSizeBPs))
                chromLength = posList[-1]
                print("chromosome length = {}".format(chromLength))
                fws=float(windowSizeBPs)
                iws=int(windowSizeBPs)
                lastEdge = int(math.ceil(chromLength / fws)) * iws
                allEdges = range(0, lastEdge, iws)
                lefts=allEdges[0:-1]
                rights=allEdges[1:]
                print("Breaking into {} windows from 0 to {}:".format(len(lefts), lastEdge))
                print(allEdges)
                for i in range(len(lefts)):
                    windowNumber='w{}'.format(i)
                    fromBp = lefts[i]
                    toBp = rights[i]
                    if raxTrees:
                        outFasta = "{}/{}_{}.fasta".format(directory, vcfPrefix, windowNumber)
                        treeCommand = "/work/02260/grovesd/stampede2/standard-RAxML/raxmlHPC-PTHREADS -s {} -n {} -m GTRCAT -f d -p 12345 -T 1".format(outFasta, outFasta.rstrip(".fasta").split("/")[1])
                        sliceCommand = "vcftools --vcf {} --chr {} --from-bp {} --to-bp {} --recode -c | vk phylo fasta - > {} && {}\n".format(vcf, chrom, fromBp, toBp, outFasta, treeCommand)
                    else:
                        sliceCommand = "vcftools --vcf {} --chr {} --from-bp {} --to-bp {} --recode -c | vk phylo fasta - > {}/{}_{}.fasta\n".format(vcf, chrom, fromBp, toBp, directory, vcfPrefix, windowNumber)
                    out.write(sliceCommand)
                    if quickTrees:
                        treeCommand = "vcftools --vcf {} --chr {} --from-bp {} --to-bp {} --recode -c | vk phylo tree upgma - > {}/{}_{}.newick\n".format(vcf, chrom, fromBp, toBp, treeDir, vcfPrefix, windowNumber)
                        out.write(treeCommand)
                    windowList.append([windowNumber, str(fromBp), str(toBp)])

            windowInfoFile="{}/{}_w{}_s{}_windowData.tsv".format(directory, chrom, windowSize, slideLength)
            print("\nWriting out the window data for chromosome {} to the file {}".format(chrom, windowInfoFile))
            with open(windowInfoFile, 'w') as out2:
                for win in windowList:
                    out2.write("\t".join(win) + "\n")




#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))        


