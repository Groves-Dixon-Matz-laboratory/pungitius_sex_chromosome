#!/usr/bin/env python
## haploidize_phased_vcf.py
##last updated 1-31-18
#Reads in a phased VCF and outputs a new VCF with a column for each chromosome
#the two columns for a sample have suffixes '_A' and '_B'


#import modules
import time
import argparse
from sys import argv
from sys import exit
import numpy as np






#Functions
def read_file(infileName, outfileName):
    '''Function to read in the VCF and split the phased data into haploid sample columns
    '''
    lineNumber = 0
    preHeadList = []
    with open(infileName, 'r') as infile:
        with open(outfileName, 'w') as out:
            for line in infile:
                lineNumber += 1
                if line[0:2] == "##":
                    out.write(line)
                    continue
                elif line[0][0] == "#":
                    line = line.strip('\n').split("\t")
                    header=line
                    sampleStart = header.index("FORMAT") + 1
                    sampleList = header[sampleStart:]
                    print("\nSamples Found in VCF:")
                    print(sampleList)
                    hapSampleList = []
                    for s in sampleList:
                        hapSampleList.append(s+"_A")
                        hapSampleList.append(s+"_B")
                    print("\nNew 'haploid' sample List:")
                    print(hapSampleList)
                    headerString = "\t".join(line[0:sampleStart] + hapSampleList) + "\n"
                    # out.write("##Haploidized by haploidize_phased_vcf.py\n")
                    out.write(headerString)
                    continue
                else:
                    line = line.strip('\n').split("\t")
                    locusData = line[0:sampleStart]
                    genotypeData = line[sampleStart:]
                    newGenoData = []
                    for g in genotypeData:
                        haps = g.split("|")
                        newGenoData.append(haps[0])
                        newGenoData.append(haps[1])
                    outString = "\t".join(locusData + newGenoData) + "\n"
                    out.write(outString)



#main
def main():
    Start_time = time.time() ##keeps track of how long the script takes to run
    

    ##Set Up Argument Parsing
    Description = '''
    Description:

    Haploidize a phased vcf (a new column for each chromosome).

    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input', help = 'The the input file vcf file phased by beagle')
    parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')
    args = parser.parse_args()


    infileName = args.input
    outfileName = args.out

    read_file(infileName, outfileName)
    
    #return time to run
    Time = time.time() - Start_time
    print('\nTime took to run: {}'.format(Time))        



if __name__ == "__main__":
   main()   


