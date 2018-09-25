#!/usr/bin/env python
##vcf_window_wrapper.py
##written 2/7/18 
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit



##############################
###### DEFINE FUNCTIONS ######
##############################
def read_hvgs():
    '''Read in the hgvs files. Keep their file names and
    the number of variants expected for each provean output file.

    Output a dictionary with the variant amino acids tested for
    each gene. These are for comparison with the output files
    to make sure everything completed.'''
    print("\nReading in HGVS files...")
    hvgsDict = {}
    stopCodonDict = {}
    geneNameList = []
    hasStopList = []
    fCount = 0
    for infileName in varFileList:
        nsplit = infileName.split("_")
        geneName = nsplit[-3] + "_" + nsplit[-2]
        with open(infileName, 'r') as infile:
            gVarList = []
            for line in infile:
                if "*" in line:
                    hasStopList.append(geneName)
                    try:
                        stopCodonDict[geneName].append(line)
                    except KeyError:
                        stopCodonDict[geneName] = [line.strip("\n")]
                else:
                    gVarList.append(line.strip("\n"))
        geneNameList.append(geneName)
        hvgsDict[geneName] = gVarList
        fCount += 1
    print("{} total files".format(fCount))
    print("{} total genes has a premature stop codon".format(len(hasStopList)))
    return(geneNameList, hvgsDict, hasStopList, stopCodonDict)



def read_provean_results():
    print("\nReading in provean results files...")
    pGeneList = []
    completedList = []
    completedGenes = []
    incompleteVarList = []
    noVarsList = []
    failedList = []
    pResDict = {}
    for infileName in infileList:
        nsplit = infileName.split("_")
        geneName = nsplit[0] + "_" + nsplit[1]
        pGeneList.append(geneName)
        with open(infileName, 'r') as infile:
            varList = []
            scoreList = []
            record = False
            wroteScores = False
            completed = False
            for line in infile:
                if "VARIATION" in line:
                    if "SCORE" in line:
                        record = True
                        wroteScores = True
                else:
                    if record:
                        line=line.strip("\n").split()
                        varList.append(line[0])
                        scoreList.append(line[1])
                    else:
                        continue
            #now have recorded set of vars and scores from provean results file
            #check against the expected number from input
            inVars = hvgsDict[geneName]
            match = inVars==varList
            # print("\n------------")
            # print(geneName)
            # print('all vars: {}'.format(match))
            
            #now assign the file to category based on completeness
            if match:
                completedList.append(infileName)
                completedGenes.append(geneName)
                pResDict[geneName] = [varList, scoreList]
            else:
                failedList.append(geneName)
                if wroteScores:
                    incompleteVarList.append(infileName)
                else:
                    noVarsList.append(infileName)
    #summarize results
    missingProveanList = []
    for g in geneNameList:
        if g not in pGeneList:
            missingProveanList.append(g)
    print("\n{} total provean output files examined".format(len(pGeneList)))
    print("{} were complete, with scores for each varaint in associated hgvs file".format(len(completedList)))
    print("{} were incomplete, with some scores, but not the full set".format(len(incompleteVarList)))
    print("{} failed to record any variant scores".format(len(noVarsList)))
    print("{} did not have and provean result file at all".format(len(missingProveanList)))
    if ( len(completedList) + len(incompleteVarList) + len(noVarsList) + len(missingProveanList) ) == len(geneNameList):
        print("Good. All files account for.")
    else:
        print("Error. Missing some files.")
    print("Writing out failed runs as failed_runs.txt")
    with open('failed_runs.txt', 'w') as out:
        for fail in missingProveanList + failedList:
            out.write(fail + "\n")
    return(completedGenes, pResDict)



def read_gtf():
    coordDict = {}
    with open(gtfInput, 'r') as infile:
        for line in infile:
            line = line.strip("\n").split("\t")
            lineChrom = line[0]
            if chrom:
                if lineChrom != chrom:
                    continue
            featureType = line[2]
            if featureType == 'gene':
                dat = line[8]
                geneId = dat.split("gene_id ")[1].split(";")[0]
                tranId = dat.split("transcript_id ")[1]
                geneName = "{}_{}".format(geneId, tranId)
                start = line[3]
                stop = line[4]
                coordDict[geneName] = [start, stop]
            #at least one gene in gtf appears to lack a 'gene' line
            #here add coordinates based on 'transcript' in case this occured
            if featureType == 'transcript':
                dat = line[8]
                geneId = dat.split("gene_id ")[1].split(";")[0]
                tranId = dat.split("transcript_id ")[1]
                geneName = "{}_{}".format(geneId, tranId)
                start = line[3]
                stop = line[4]
                try:
                    coordDict[geneName]
                except KeyError:
                    coordDict[geneName] = [start, stop]

    return coordDict




def output_results():
    varsWritten = 0
    stopsWritten = 0
    with open(outfileName, 'w') as out:
        out.write("gene\tstart\tstop\tvar\tscore\tstopCodon")
        for g in completedGenes:
            vdat = pResDict[g]
            varList = vdat[0]
            scoreList = vdat[1]
            cdat = coordDict[g]
            s = cdat[0]
            e = cdat[1]
            for i in range(len(varList)):
                v=varList[i]
                score=scoreList[i]
                outList = [g, s, e, v, score, "FALSE"]
                outString = "\t".join(outList)
                out.write("\n{}".format(outString))
                varsWritten += 1
            #now add stop codons if any
            try:
                stopDat = stopCodonDict[g]
                for v in stopDat:
                    v=v.strip("\n")
                    outList = [g, s, e, v, "NA", "TRUE"]
                    outString = "\t".join(outList)
                    out.write("\n{}".format(outString))
                    stopsWritten += 1
            except KeyError:
                continue
    print("\n{} total variants written".format(varsWritten))
    print("{} total stops written".format(stopsWritten))
    print("results written to: {}".format(outfileName))










##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'provean_res_files', nargs="+", help = 'Glob to the provean files')
    parser.add_argument('-hgvs', required = True, dest = 'the_variant_files', nargs="+", help = 'Glob to the hgvs files used to run provean')
    parser.add_argument('-gtf', required = True, dest = 'input_gtf', help = 'the input gtf file. Used for gathering gene coordinates.')
    parser.add_argument('-chr', required = False, default=False, dest = 'chromosome', help = 'optional target chromosome. All others ignored')
    parser.add_argument('-o', required = True, dest = 'output_file', help = 'output name')


    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infileList = args.provean_res_files
    varFileList = args.the_variant_files
    outfileName = args.output_file
    gtfInput = args.input_gtf
    chrom = args.chromosome

    #--- RUN FUNCTIONS ---#
    coordDict = read_gtf()
    geneNameList, hvgsDict, hasStopList, stopCodonDict = read_hvgs()
    completedGenes, pResDict = read_provean_results()
    output_results()
    print("\nDone.\n")


