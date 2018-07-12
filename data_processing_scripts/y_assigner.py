#!/usr/bin/env python
## y_assigner.py
##last updated 4-16-18


#import modules
import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio import Phylo
import sys

debug = False


#Functions
def read_trees(infileName, monoGroupList, monoSampleNameList, outfileName, windowList, rooted):
    '''Function to read in a file as a list of lists
    '''
    with open(labeledTreesOut, 'w') as calledYOutFile:
        print("\nAssuming an unrooted tree. (will try to detect monopyly in both directions)")
        expectedForSex = len(monoGroupList) / 2
        #iterate through the trees in the newick file
        treeCount = 0
        treeIndex = -1
        treeRecordList = []
        biggestMonoList = []
        switchesNeeded = []
        windowRecordList = []
        SDRboolList = [] ##keep track of whether tree passes rule: at least one chromosome from the male sample set is found in a monophyletic group of males
        trees = Phylo.parse(infileName, 'newick')
        totTrees = 0
        for t in trees: totTrees += 1
        trees = Phylo.parse(infileName, 'newick')
        print("\nTotal number of trees in tree file = {}".format(totTrees))
        if totTrees != len(windowList):
            print("ERROR. The number of trees does not match the window list")
        else:
            print("Tree Number matches expectation based on window file.")
        for tree in trees:
            treeCount+=1
            treeIndex+=1
            treeWindow = windowList[treeIndex]
            biggestMonoSize = 0
            atLeastOneChrom = False
            biggestYlist = []
            smallestSDRlike = []
            nodeList = tree.get_nonterminals()
            nodeCount = 0
            for n in nodeList: nodeCount += 1
            if debug:
                print("-------------------------")
                print("Tree #{}:".format(treeCount))
                print("Window = {}".format(treeWindow))
                print("{} terminal taxa".format(tree.count_terminals()))
                print("{} internal nodes".format(nodeCount))
            terminalList = [term.name for term in tree.get_terminals()]
            #make sure the taxa given in input file are in the tree
            if treeCount % 10 == 0:
                print("Reading tree number {}...".format(treeCount))
            check_tree_has_group(monoGroupList, terminalList)



            nonTerminals = tree.get_nonterminals()
            for nt in nonTerminals:
                isMono = True
                monoSize = 0
                nodeTerminals = nt.get_terminals()
                
                #check the internal node's leaves to see if monophyletic for specific group
                leafList = [t.name for t in nodeTerminals]
                biggestMonoSize, atLeastOneChrom, biggestYlist, smallestSDRlike = check_leaf_set(leafList, monoGroupList, monoSampleNameList, biggestMonoSize, atLeastOneChrom, biggestYlist, smallestSDRlike)



                #if the tree is unrooted, you need to check the opposite direction
                if not rooted:
                    antiLeafList = []
                    for treeLeaf in terminalList:
                        if treeLeaf not in leafList:
                            antiLeafList.append(treeLeaf)
                    biggestMonoSize, atLeastOneChrom, biggestYlist, smallestSDRlike = check_leaf_set(leafList, monoGroupList, monoSampleNameList, biggestMonoSize, atLeastOneChrom, biggestYlist, smallestSDRlike)
            treeRecordList.append(str(treeCount))
            biggestMonoList.append(str(biggestMonoSize))
            windowRecordList.append(windowList[treeIndex])
            switchesNeeded.append(str(expectedForSex - biggestMonoSize))
            SDRboolList.append(str(atLeastOneChrom))

            #append results for this tree to the output lists
            if debug:
                print("Biggest monophyletic group size = {}".format(biggestMonoSize))
                print("Expcected for a sex chromosome = {}".format(expectedForSex))
                print("Total switches needed = {}".format(expectedForSex - biggestMonoSize))
                print("Monophyletic clade with at least one chrom from each male? ('SDR-like')  = {}".format(atLeastOneChrom))
                print("Smallest SDR-like clade size = {}".format(smallestSDRlike))

            #Add yness scores to all the positions in this tree
            calledYs = score_Yness(treeWindow, biggestYlist, smallestSDRlike, atLeastOneChrom)
            newTree = tree
            nonTerminals = newTree.get_terminals()
            for nt in nonTerminals:
                if nt.name in monoGroupList:
                    nt.name = nt.name + "_M"
                if nt.name.rstrip("_M") in calledYs:
                    nt.name = nt.name + "__Y"
            # print('\n\n\nCalled Ys and Tree for Window#{} with labeles:'.format(treeWindow))
            # print(calledYs)
            if perfectOnly:
                if len(calledYs) == len(monoSampleNameList):
                    Phylo.write(newTree, calledYOutFile, "newick")
            else:
                Phylo.write(newTree, calledYOutFile, "newick")
            # print("\n\n")




    #output the results
    # output(outfileName, treeRecordList, windowRecordList, biggestMonoList, expectedForSex, switchesNeeded, SDRboolList)
    return treeRecordList, windowRecordList, biggestMonoList, switchesNeeded, atLeastOneChrom



def score_Yness(treeWindow, biggestYlist, smallestSDRlike, atLeastOneChrom):
    # print("Scoring Yness for tree made from window {}".format(treeWindow[0]))
    wLeft = int(treeWindow[1])
    wRight = int(treeWindow[2])
    posArr = np.array(posList)
    wArr = [pos for pos in posArr if pos >= wLeft and pos <= wRight]
    if atLeastOneChrom:
        # print("This tree had an SDR-like clade")
        # print("Smallest SDR-like clade:")

        #for each sample add which chromosome is Y-like to each position wihtin this window
        for maleChrom in smallestSDRlike:
            sample = maleChrom.split("_")[0]
            for pos in wArr:
                ydict[pos][sample].append(maleChrom)
        toReturn = smallestSDRlike
    else:
        # print("This tree did NOT have and SDR-like clade")
        # print("Biggest monophyletic group:")
        if len(biggestYlist) > monoCutoff*len(monoSampleNameList):
            for maleChrom in biggestYlist:
                sample = maleChrom.split("_")[0]
                for pos in wArr:
                    ydict[pos][sample].append(maleChrom)
            toReturn = biggestYlist
        else:
            toReturn = 'none called'
    return toReturn



def check_leaf_set(aLeafList, monoGroupList, monoSampleNameList, biggestMonoSize, atLeastOneChrom, biggestYlist, smallestSDRlike):
    """Funciton to compare the leaves of an interal node against the monogroup list to check if it is monophyletic for them.
    Also checks if the clade is 'SDR-like', meaning that at least one chromosome from each sample is present.
    If it is monophyletic, keep the largest clade as biggestYlist
    If it is 'SDR-like', keep the smallest clade.
    The smallest 'SDR-like' clade is assumed to be the most specific monophyletic clade that still has
    each male represented. Ideally the size will be equal to the number of male samples.
    (ie try to keep our mostly X chromosomes with contaminating Y alleles)
    """
    monoSize=0
    isMono=True
    for l in aLeafList:
        if l not in monoGroupList:
            isMono = False
            break
        else:
            monoSize += 1
    if isMono:
        if monoSize > biggestMonoSize:
            biggestMonoSize = monoSize
            biggestYlist = [leaf for leaf in aLeafList] #record the leaves for biggest monophyletic
        selectSampleList = [x[0:-2] for x in aLeafList]
        # print("Found a monophyletic group")
        # print("Samples in this leaf list:")
        # print(selectSampleList)
        sInBoth = []
        for s in monoSampleNameList:
            if s in selectSampleList:
                sInBoth.append(s)
        if len(sInBoth) == len(monoSampleNameList):
            atLeastOneChrom = True
            # print("This clade counts as an SDR")
            newLeafList = [leaf for leaf in aLeafList]
            if len(newLeafList) > len(smallestSDRlike):
                smallestSDRlike = newLeafList
    return biggestMonoSize, atLeastOneChrom, biggestYlist, smallestSDRlike



def read_monogroup(monoGroupFile):
    """Funciton to read the monoGroupFile into a list"""
    print("\nReading in group file to test for biggest monophyletic group...")
    monoGroupList = []
    monoSampleNameList = []
    with open(monoGroupFile, 'r') as infile:
        for line in infile:
            sampleName = line.strip("\n")
            monoSampleNameList.append(sampleName)
            monoGroupList.append(sampleName + "_A")
            monoGroupList.append(sampleName + "_B")
    print("Found {} samples in monoGroupFile".format(len(monoGroupList)))
    print("Will look for the following terminal taxa in the trees:")
    for i in monoGroupList: print i
    print("If these fail to appear in a tree an error will be thrown")
    return monoGroupList, monoSampleNameList


def check_tree_has_group(monoGroupList, terminalList):
    """Funciton to check that the individuals in the monoGroupFile
    can be found in a tree. If not returns the missing individuals."""
    missingList = []
    for t in monoGroupList:
        if t not in terminalList:
            missingList.append(t)
    if len(missingList)>0:
        print("Error. Following taxa from monoGroupFile are missing from this tree:")
        for x in missingList: print x
        exit("Please check trees and input file. Quitting.")
    else:
        # print("All monoGroupFile individuals found in tree")
        pass


def read_window(windowFile):
    windowList = []
    with open(windowFile, 'r') as infile:
        for line in infile:
            windowList.append(line.strip("\n").split())
    return windowList

def output(outfileName, treeRecordList, windowRecordList, biggestMonoList, expectedForSex, switchesNeeded, SDRboolList):
    """Funciton to output the results to a file"""
    print("\nOutputting results to file {}...".format(outfileName))
    if len(treeRecordList) != len(biggestMonoList) or len(switchesNeeded) != len(treeRecordList) != len(SDRboolList):
        exit("Error! Bug in script.")
    with open(outfileName, 'w') as out:
        header = ['treeNumber', 'windowNumber', 'BIN_START', 'BIN_END', 'biggestMonoPhyletic', 'expected', 'switchesNeeded', 'SDRbool']
        out.write("\t".join(header))
        for i in range(len(treeRecordList)):
            outList = [treeRecordList[i], windowRecordList[i], biggestMonoList[i], str(expectedForSex), switchesNeeded[i], SDRboolList[i]]
            out.write("\n" + "\t".join(outList))


def read_vcf(vcfFile):
    """read in the vcf file"""
    print("Reading VCF file {}...".format(vcfFile))
    posList = []  #keep track of all the positions in the VCF
    with open(vcfFile, 'r') as infile:
        for line in infile:
            #skip header lines
            if line[0] == "#":
                continue
            #count up variants
            line=line.split()
            pos=int(line[1])
            posList.append(pos)
    print("\nFound {} variants in the VCF file {}".format(len(posList), vcfFile))
    return(posList)


def setup_y_score_dict(posList, monoSampleNameList):
    """Set up dictionary to keep score of which chromosome (A or B) from each sample
    in the monogroupList is called as the Y in each window. If a sliding window
    was used, each site will get multiple calls."""
    ydict = {}
    for p in posList:
        ydict[p] = {}
        for male in monoSampleNameList:
            ydict[p][male] = []
    return(ydict)




def finalize_yassignment():
    print("\nFinalizing results...")
    with open(outYFile, 'w') as out:
        out.write("\t".join(["position", "sample", "AisY", "BisY"]))
        for p in posList:
            for m in monoSampleNameList:
                letters = [x.split("_")[1] for x in ydict[p][m]]
                aCount = letters.count("A")
                bCount = letters.count("B")
                if aCount + bCount != len(letters):
                    exit("Error! Counts in ydict didn't make sense")
                outString = "\n{}\t{}\t{}\t{}".format(p, m, aCount, bCount)
                out.write(outString)
    if perfectOnly:
        print("\nTrees perfect Y assignments saved as {}".format(labeledTreesOut))
    else:
        print("\nAll trees with Y chomosomes labeled saved as {}".format(labeledTreesOut))








#main
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
Description = '''
Description:

'''

parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-trees', required = False, dest = 'input_trees', help = 'The the input file')
parser.add_argument('-rooted', required = False, default=False, dest = 'rooted', help = 'Change to True if tree is rooted')
parser.add_argument('-i', required = True, dest = 'input_group_file', help = 'File of the sample names to look for a monophyletic group of (without _A and _B appended)')
parser.add_argument('-w', required = True, dest = 'windows_file', help = 'Name for the window annotations for these trees')
parser.add_argument('-vcf', required = True, dest = 'vcf_file', help = 'The VCF file used to generate the windows/trees')
parser.add_argument('-o', required = False, dest = 'output_name', help = 'Name for output file')
parser.add_argument('-yout', required = True, dest = 'y_output_name', help = 'Name for output file giving scores for which chromosome is Y at each position')
args = parser.parse_args()


infileName = args.input_trees
rooted = args.rooted
monoGroupFile = args.input_group_file
windowFile = args.windows_file
outfileName = args.output_name
vcfFile = args.vcf_file
outYFile = args.y_output_name

labeledTreesOut = infileName.rstrip(".newick") + "_Ylabeled.newick"
monoCutoff = 0.75
perfectOnly = True


#RUN FUNCTIONS
posList = read_vcf(vcfFile)
windowList = read_window(windowFile)
monoGroupList, monoSampleNameList = read_monogroup(monoGroupFile)
ydict = setup_y_score_dict(posList, monoSampleNameList)
read_trees(infileName, monoGroupList, monoSampleNameList, outfileName, windowList, rooted)
finalize_yassignment()



#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))        



