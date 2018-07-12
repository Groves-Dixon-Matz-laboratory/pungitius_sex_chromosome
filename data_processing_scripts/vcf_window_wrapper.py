#!/usr/bin/env python
##vcf_window_wrapper.py
##written 2/7/18 
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
import numpy as np
import subprocess
import math
import os
from multiprocessing import Pool
import pandas as pd



##############################
###### DEFINE FUNCTIONS ######
##############################

def read_sex(sexTable):
    with open(sexTable, 'r') as infile:
        sexList = [line.strip() for line in infile if line.rstrip()]
    print(sexList)
    return(sexList)

def call_vcftools_fst():
    """Call VCFtools to run Fst on windows."""
    print("\nRunning Windowed Fst...")
    stdO=open("{}/windowFst.o".format(logDir), 'w')
    stdE=open("{}/windowFst.e".format(logDir), 'w')
    fstOut = "{}/{}".format(resultsDir, prefix)
    cmd = ['vcftools', "--vcf", inputVcf, "--weir-fst-pop", maleFile, "--weir-fst-pop", femaleFile, "--fst-window-size", str(windowSize), "--fst-window-step", str(slideSize), "--out", fstOut]
    subprocess.Popen(cmd, stdout=stdO, stderr=stdE)
    stdO.close()
    stdE.close()
    print("Done.")

def call_vcftools_split_sex():
    print("\nstarting splitting VCF by sex...")
    stdO=open("{}/splitSex.o".format(logDir), 'w')
    stdE=open("{}/splitSex.e".format(logDir), 'w')
    maleOut = "{}/{}_male".format(winDir, prefix)
    femaleOut = "{}/{}_female".format(winDir, prefix)
    mcmd = ['vcftools', "--vcf", inputVcf, "--keep", maleFile, "--recode", "--out", maleOut]
    fcmd = ['vcftools', "--vcf", inputVcf, "--keep", femaleFile, "--recode", "--out", femaleOut]
    splitMaleProc = subprocess.Popen(mcmd, stdout=stdO, stderr=stdE)
    splitFemaleProc = subprocess.Popen(fcmd, stdout=stdO, stderr=stdE)
    maleVcf = maleOut + ".recode.vcf"
    femaleVcf = femaleOut + ".recode.vcf"
    return splitMaleProc, maleVcf, splitFemaleProc, femaleVcf, stdE, stdO

def call_vcftools_split_sex_for_snp_density():
    print("\nstarting splitting VCF by sex...")
    stdO=open("{}/splitSexForSnp.o".format(logDir), 'w')
    stdE=open("{}/splitSexForSnp.e".format(logDir), 'w')
    maleOut = "{}/{}_maleForSnp".format(winDir, prefix)
    femaleOut = "{}/{}_femaleForSnp".format(winDir, prefix)
    mcmd = ['vcftools', "--vcf", inputVcf, "--keep", maleFile, "--maf", "0.1", "--recode", "--out", maleOut]
    fcmd = ['vcftools', "--vcf", inputVcf, "--keep", femaleFile, "--maf", "0.1", "--recode", "--out", femaleOut]
    splitMaleProc = subprocess.Popen(mcmd, stdout=stdO, stderr=stdE)
    splitFemaleProc = subprocess.Popen(fcmd, stdout=stdO, stderr=stdE)
    maleVcf = maleOut + ".recode.vcf"
    femaleVcf = femaleOut + ".recode.vcf"
    return splitMaleProc, maleVcf, splitFemaleProc, femaleVcf, stdE, stdO


def call_vcftools(vcfFile, vcfToolsArgumentList, callTag):
    print("\nRunning {} for VCF file {}...".format(callTag, vcfFile))
    stdO=open("{}/{}.o".format(logDir, callTag), 'w')
    stdE=open("{}/{}.e".format(logDir, callTag), 'w')
    subPrefix = vcfFile.split("/")[-1].rstrip(".vcf")
    subOut = "{}/{}".format(resultsDir, subPrefix)
    # print("writing results to {}".format(subOut))
    cmd = ['vcftools', "--vcf", vcfFile, "--out", subOut] + vcfToolsArgumentList
    p=subprocess.Popen(cmd, stdout=stdO, stderr=stdE)
    return(p, stdO, stdE)
    

def format_af_df(frqFile):
    """reformat an allele frequency dataframe output from VCFtools --freq"""
    dat = pd.read_csv(frqFile, sep="\t", index_col=False, header = 2, names=['chr', 'pos', 'nAllele', 'nChr', 'af1', 'af2'])
    af1=dat["af1"]
    af2=dat["af2"]
    dat["pos"] = pd.to_numeric(dat["pos"], downcast="integer")
    dat["a1"] = af1.str.split(":").str.get(0)
    dat["a2"] = af2.str.split(":").str.get(0)

    dat["f1"] = pd.to_numeric(af1.str.split(":").str.get(1), errors = 'coerce')
    dat["f2"] = pd.to_numeric(af2.str.split(":").str.get(1), errors = 'coerce')
    return(dat)


def getr2(sub1, sub2):
    ip1 = sub1['f1']
    ip2 = sub2['f1']
    isPoly = (ip1 >= 0.1) & (ip1 <= 0.9)
    p1 = ip1.ix[isPoly]
    p2 = ip2.ix[isPoly]
    r2=(p1-p2)**2 / (p1*(1-p1))
    return(r2)


def allele_frequency_stats(femaleVcf, maleVcf, windowList):
    """Function to gather statistics in windows that VCF doesn't do automatically.
    Reads in the allele frequency output files for males and females.
    Also reads in the --site-depth output files for the males and females."""
    precenceCutoff = 0.4
    print("\nGetting stats from allele frequency data...")

    #SET UP THE INPUT FILES
    #assign prefixes
    femalePrefix = femaleVcf.split("/")[-1].rstrip(".vcf")
    malePrefix = maleVcf.split("/")[-1].rstrip(".vcf")
    #assign allele freq filename
    femaleFrqFile = "{}/{}.frq".format(resultsDir, femalePrefix)
    maleFrqFile = "{}/{}.frq".format(resultsDir, malePrefix)
    #assign depth filename
    femaleDepthFile = "{}/{}.ldepth".format(resultsDir, femalePrefix)
    maleDepthFile = "{}/{}.ldepth".format(resultsDir, malePrefix)
    #read in allele freq
    fdat=format_af_df(femaleFrqFile)
    mdat=format_af_df(maleFrqFile)
    #read in depth
    fDdat = pd.read_csv(femaleDepthFile, sep="\t")
    mDdat = pd.read_csv(maleDepthFile, sep="\t")
    fDdat["pos"] = fDdat["POS"] #just to match 
    mDdat["pos"] = mDdat["POS"]

    #Double-check that the files agree as they should
    print("\nDouble-check the positions in the allele frequency files match up:")
    print(sum(fdat["pos"]==mdat["pos"])==mdat.shape[0])
    if not (sum(fdat["pos"]==mdat["pos"])==mdat.shape[0]):
        exit("\nError. The allele frequency tables do not match between sexes. There is a bug somewhere.")
    if not (sum(fDdat["POS"]==mDdat["POS"])==fDdat.shape[0]):
        exit("\nError. The site depth tables do not match between sexes. There is a bug somewhere.")



    #get stats for windows
    dxyList = []
    mDepthList = []
    fDepthList = []
    propFemaleSpecificList = []
    propMaleSpecificList = []
    propEitherSpecificList = []
    meanDiffA1List = []
    meanDiffA2List = []
    fr2List = []
    mr2List = []
    leftList = []
    rightList = []
    for i in range(len(windowList)):
        l=int(windowList[i][0])
        r=int(windowList[i][1])
        leftList.append(l)
        rightList.append(r)
        msub = mdat.ix[(mdat["pos"] >=l) & (mdat["pos"] <= r)]
        fsub = fdat.ix[(fdat["pos"] >=l) & (fdat["pos"] <= r)]
        mDsub = mDdat.ix[(mDdat["pos"] >=l) & (mDdat["pos"] <= r)]
        fDsub = fDdat.ix[(fDdat["pos"] >=l) & (fDdat["pos"] <= r)]
        mtot=msub.shape[0]
        ftot=fsub.shape[0]
        if mtot != ftot:
            exit("\nError, windows are not same between male and female .frq files.\nThere is a bug somewhere.\n")
        
        #handle possibility that this window does not have any SNPs in it
        if mtot == 0:
            dxyList.append('no_snps_in_this_window')
            propFemaleSpecificList.append('no_snps_in_this_window')
            propMaleSpecificList.append('no_snps_in_this_window')
            propEitherSpecificList.append('no_snps_in_this_window')
            meanDiffA1List.append('no_snps_in_this_window')
            meanDiffA2List.append('no_snps_in_this_window')
            fr2List.append('no_snps_in_this_window')
            mr2List.append('no_snps_in_this_window')
            continue

        
        #GET PERCENT SEX-SPECIFIC LOCI PER WINDOW
        #get bool arrays for frequency == 0 in each sex:
        faf1_0 = fsub["f1"] == 0  #female freqeuncy for allele 1 is zero
        faf2_0 = fsub["f2"] == 0  #female freqeuncy for allele 2 is zero
        maf1_0 = msub["f1"] == 0  #male frequency for allele 1 is zero
        maf2_0 = msub["f2"] == 0  #male frequency for allele 2 is zero

        #get bool arrays for precent of allele in each sex based on cutoff
        faf1_p = fsub["f1"] >= precenceCutoff  #allele 1 present in females
        faf2_p = fsub["f2"] >= precenceCutoff  #allele 2 present in females
        maf1_p = msub["f1"] >= precenceCutoff  #allele 1 present in males
        maf2_p = msub["f2"] >= precenceCutoff  #allele 2 present in males

        #make arrays for alleles found in one sex but not the other
        fSpecific1 = faf1_p & maf1_0  #allele 1 present in females but absent in males
        fSpecific2 = faf2_p & maf2_0  #allele 2 present in females but absent in males
        mSpecific1 = maf1_p & faf1_0  #allele 1 present in males but absent in females
        mSpecific2 = maf2_p & faf2_0  #allele 2 present in males but absent in females


        #make either sex-specific
        eitherSpecific = fSpecific1 | fSpecific2 | mSpecific1 | mSpecific2

        #get proportions
        propFemaleSpecific = float(sum(fSpecific1) + sum(fSpecific2)) / float(ftot)
        propMaleSpecific = float(sum(mSpecific1) + sum(mSpecific2)) / float(mtot)

        propEitherSpecific = float(sum(eitherSpecific)) / ( ( float(ftot) + float(mtot) ) / 2 )

        propFemaleSpecificList.append(propFemaleSpecific)
        propMaleSpecificList.append(propMaleSpecific)
        propEitherSpecificList.append(propEitherSpecific)




        #GET MEAN MALE AND FEMALE SITE-DEPTH SUM
        #REMOVED THIS
        # mnMaleD = np.mean(mDsub["SUM_DEPTH"])
        # mnFemaleD = np.mean(fDsub["SUM_DEPTH"])
        # mDepthList.append(mnMaleD)
        # fDepthList.append(mnFemaleD)


        #GET DXY
        #dxy = (P1)(1-p2) + (p2)(1-p1) / #bp
        p1 = msub["f1"]
        p2 = fsub["f1"]
        q1 = 1.0 - p1
        q2 = 1.0 - p2
        diffs = p1*q2 + p2*q1
        dxy = float(sum( diffs )) / float(windowSize)
        dxyList.append(dxy)


        #GET MEAN ALLELE FREQUENCY DIFFERENCE
        diffA1 = abs(msub["f1"] - fsub["f1"])
        diffA2 = abs(msub["f2"] - fsub["f2"])
        meanDiffA1List.append(np.mean(diffA1))
        meanDiffA2List.append(np.mean(diffA2))

        #GET r2 SEX-LINKAGE STATISTIC FOR MALE AND FEMALE
        fr2 = getr2(fsub, msub)
        mr2 = getr2(msub, fsub)
        fr2List.append(np.mean(fr2))
        mr2List.append(np.mean(mr2))
    chromList = [chrom]*len(leftList)
    resultListList = [chromList, leftList, rightList, dxyList, propFemaleSpecificList, propMaleSpecificList, propEitherSpecificList, meanDiffA1List, meanDiffA2List, fr2List, mr2List]
    colNames = ['CHROM', 'BIN_START', 'BIN_END', 'Dxy', 'pFemaleSpecific', 'pMaleSpecific', 'pSexSpecific', 'meanDiffA1', 'meanDiffA2', 'female_r2', 'male_r2']
    print("Results")
    for ci in range(len(colNames)):
        print("{} : {}".format(colNames[ci], len(resultListList[ci])))
    res = pd.DataFrame(
        {'CHROM': chromList,
        'BIN_START': leftList, 
        'BIN_END': rightList,
        'Dxy': dxyList,
        'pFemaleSpecific': propFemaleSpecificList, 
        'pMaleSpecific': propMaleSpecificList,
        'pSexSpecific': propEitherSpecificList,
        'meanDiffA1': meanDiffA1List, 
        'meanDiffA2': meanDiffA2List, 
        'female_r2': fr2List, 
        'male_r2': mr2List})[colNames]
    print("------------------------")
    print("\nDone calculating allele frequency stats")
    print("results:")
    print(res.head())
    print(res.shape)
    res.to_csv(windowFreqResultsFile, sep="\t")






def window_vcf(inputVcf, windowSize, slideSize, chrom, chromLength):
    '''Take the window size, slide size, and chromosome length to build the lists
    of left and right window boundaries. Output is a list with length equal to the
    number of windows, where each element is a sublist with a paired of left and right
    coordinates.
    eg for window size = 1000 and slide = 1000: windowList = [[1, 1000], [1001, 2000]...]
    '''
    print("Setting up windows {}bp wide with slide length = {}...".format(inputVcf, windowSize, slideSize))
    lastEdge = int(math.ceil(float(chromLength)/float(windowSize))*windowSize)
    lefts = range(1, lastEdge, slideSize)
    rights = range(windowSize, lastEdge, slideSize)
    rights.append(lastEdge)
    windowList = [[lefts[i], rights[i]] for i in range(len(lefts))]
    print("\nmaking windows for vcf {}".format(inputVcf))
    print("{} Total windows from {} to {}".format(len(lefts), lefts[0], rights[-1]))
    prefix = inputVcf.rstrip(".vcf")
    winDir = "{}_w{}_s{}".format(prefix, windowSize, slideSize)
    if not os.path.exists(winDir):
        os.makedirs(winDir)
    return(windowList)

def call_vcftools_window(i):
    """Call VCF tools to split the vcf into a window based on left and right
    coordinates given in in the windowList at index i.
    Individual VCFs for each window are saved in the directory 'winDir' """
    stdE = open("windowSplittingLog.e", 'w')
    l=str(windowList[i][0])
    r=str(windowList[i][1])
    cmd = ['vcftools', "--vcf", inputVcf, "--chr", chrom, "--from-bp", l, "--to-bp", r, "--recode", "--out", "{}/{}_w{}_{}-{}".format(winDir, prefix, i, l, r)]
    subprocess.Popen(cmd, stdout=None, stderr=stdE)
    stdE.close()





def check_wait(process, nameTag):
    print("\nChecking if process '{}' is finished...".format(nameTag))
    if process.poll() == None:
        print("   process not finished yet, waiting...")
        process.wait()
        # while process.poll() == None:
        #     print(".")
        #     time.sleep(2)
    print("Process '{}' is complete.".format(nameTag))



def add_collumns(df1, df2, colList1, colList2, colsToAdd, namesToAdd):
    """Check that dataframes are equal for lists of columns"""
    print("\nAdding the following data columns".format(df2))
    for c in colsToAdd: print(c)
    print("Do windows from the two dataframes match?")
    allMatch = True
    failedCols = []
    tot1 = df1.shape[0]
    tot2 = df2.shape[0]
    x=df1.merge(df2, 'left', on=['BIN_START'], sort = True)
    # print("---------------------")
    # print('df2:')
    # print(df2.head())
    # print("---------------------")
    # print('m:')
    # print(x.head())
    # print("---------------------")
    # print(df1.head())
    # print('MATCH?')
    # print(df2.shape)
    # print(x.shape)
    # print(df1.shape)
    print(sum(x['BIN_START'] == df1['BIN_START']))
    print(sum(x['BIN_START'] == df1['BIN_START']) == df1.shape[0])
    for i in range(len(colsToAdd)):
        df1[namesToAdd[i]] = x[colsToAdd[i]]


def reformat_bins(df, sBin):
    df["BIN_END"] = df[sBin] + windowSize
    df[sBin] = df[sBin] + 1
    df["BIN_START"] = df[sBin]
    return(df)




def assemble_results():
    """Assemble all the results into single table"""
    print("\n\nAssembling all results into single file...")
    prefix=inputVcf.rstrip(".vcf")
    femalePrefix = prefix + "_female.recode"
    malePrefix = prefix + "_male.recode"
    piFile = "{}/{}".format(resultsDir, prefix + ".windowed.pi")

    #start data assembly with the PI data for both sexes
    dat = pd.read_csv(piFile, sep="\t")
    dat.columns = ['CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'PI_all']

    #ADD PI DATA FOR SEXES INDIVIDUALLY
    femalePiFile = "{}/{}".format(resultsDir, femalePrefix + ".windowed.pi")
    malePiFile = "{}/{}".format(resultsDir, malePrefix + ".windowed.pi")
    #female
    print("\nAdding data from file {}...".format(femalePiFile))
    dat1 = pd.read_csv(femalePiFile, sep="\t")
    print(dat.head())
    print(dat1.head())
    print(dat1["PI"].head())
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], ['PI'], ['PI_female'])
    print(dat.head())
    #male
    print("\nAdding data from file {}...".format(malePiFile))
    dat1 = pd.read_csv(malePiFile, sep="\t")
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], ['PI'], ['PI_male'])
    print(dat.head())


    #ADD FST DATA
    fstFile = "{}/{}".format(resultsDir, prefix + ".windowed.weir.fst")
    #all
    print("\nAdding data from file {}...".format(fstFile))
    dat1 = pd.read_csv(fstFile, sep="\t")
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], ['MEAN_FST'], ['MEAN_FST'])


    #ADD SNP DENSITY
    snpFile = "{}/{}".format(resultsDir, prefix + ".snpden")
    femaleSnpFile = "{}/{}".format(resultsDir, prefix + "_femaleForSnp.recode.snpden")
    maleSnpFile = "{}/{}".format(resultsDir, prefix + "_maleForSnp.recode.snpden")
    #all
    print("\nAdding data from file {}...".format(snpFile))
    dat1 = reformat_bins(pd.read_csv(snpFile, sep="\t"), "BIN_START")
    print(dat1.head())
    print(dat1.tail())
    print(dat.shape)
    print(dat1.shape)
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], ['SNP_COUNT'], ['SNP_COUNT_all'])
    #female
    print("\nAdding data from file {}...".format(femaleSnpFile))
    dat1 = reformat_bins(pd.read_csv(femaleSnpFile, sep="\t"), "BIN_START")
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], ['SNP_COUNT'], ['SNP_COUNT_female'])
    #male
    print("\nAdding data from file {}...".format(maleSnpFile))
    dat1 = reformat_bins(pd.read_csv(maleSnpFile, sep="\t"), "BIN_START")
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], ['SNP_COUNT'], ['SNP_COUNT_male'])
    print(dat.head())

    #ADD ALLELE FREQUENCY STATS
    print("\nAdding in Allele Frequency Stats...")
    dat1=pd.read_csv(windowFreqResultsFile, sep="\t", index_col=0)
    print(dat1.head())
    columnsToAdd = ['Dxy', 'pFemaleSpecific', 'pMaleSpecific', 'pSexSpecific', 'meanDiffA1', 'meanDiffA1', 'female_r2', 'male_r2']
    add_collumns(dat, dat1, ['CHROM', 'BIN_START', 'BIN_END'], ['CHROM', 'BIN_START', 'BIN_END'], columnsToAdd, columnsToAdd)
    print(dat.head())
    dat.to_csv(finalOutput, sep = "\t", na_rep = "NA")








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
    parser.add_argument('-males', required = True, dest = 'male_list', help = 'a single column text file with male sample names as they appear in the vcf')
    parser.add_argument('-females', required = True, dest = 'female_list', help = 'a single column text file with female sample names as they appear in the vcf')
    parser.add_argument('-i', required = True, dest = 'input_vcf', help = 'The the input VCF file')
    parser.add_argument('-w', required = True, dest = 'window_size', help = 'Size of the window in basepairs')
    parser.add_argument('-s', required = False, dest = 'slide_size', default = False, help = 'Distance to slide windows. Default is equal to window size')
    parser.add_argument('-l', required = True, dest = 'chrom_length', help = 'Length of the chromosome for the vcf')
    parser.add_argument('-chr', required = True, dest = 'chromosome', help = 'name of chromosome as it appears in the vcf')
    parser.add_argument('-Ncore', required = False, dest = 'N_cores', default = 1, help = 'set number of processes to run at once')
    # parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')


    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    femaleFile = args.female_list
    maleFile = args.male_list
    inputVcf = args.input_vcf
    windowSize = int(args.window_size)
    slideSize = args.slide_size
    Ncore = int(args.N_cores)
    if slideSize:
        slideSize = int(slideSize)
    else:
        slideSize = int(windowSize)
    chrom = args.chromosome
    chromLength = int(args.chrom_length)


    #--- SET UP GLOBAL VARIABLES ---#
    prefix = inputVcf.rstrip(".vcf")  #for outputs
    winDir = "VCFWRAP_SUBVCFS_{}_w{}_s{}".format(prefix, windowSize, slideSize) #for the window split vcfs
    resultsDir = "VCFWRAP_RESULTS_{}_w{}_s{}".format(prefix, windowSize, slideSize) #for results files
    logDir = "VCFWRAP_LOGFILES_{}_w{}_s{}".format(prefix, windowSize, slideSize) #for results files
    if not os.path.exists(winDir):
        os.makedirs(winDir)
    if not os.path.exists(resultsDir):
        os.makedirs(resultsDir)
    if not os.path.exists(logDir):
        os.makedirs(logDir)
    windowFreqResultsFile = resultsDir + "/" + prefix + '_alleleFequencyResults.tsv'
    finalOutput = "{}_w{}_s{}_RESULTS.tsv".format(prefix, windowSize, slideSize)


    #--- RUN FUNCTIONS ---#

    #assign the males and females
    femaleList = read_sex(femaleFile)
    maleList = read_sex(maleFile)
    nMale = len(maleList)
    nFemale = len(femaleList)
    print("\nNumber of male samples = {}".format(len(maleList)))
    print("\nNumber of female samples = {}".format(len(femaleList)))

    # assemble_results()
    # exit()


    # BEGIN SPLITTING VCF BY SEX
    splitMaleProc, maleVcf, splitFemaleProc, femaleVcf, splitE, splitO = call_vcftools_split_sex()

    ### RUN STATS ON FULL VCF
    print("\n------------------------")
    print("Running vcftools to get results from full VCF...")
    fstP, fstO, fstE = call_vcftools(inputVcf, ["--weir-fst-pop", maleFile, "--weir-fst-pop", femaleFile, "--fst-window-size", str(windowSize), "--fst-window-step", str(slideSize)], "WINDOW_FST")
    piP, piO, piE = call_vcftools(inputVcf, ["--window-pi", str(windowSize), "--window-pi-step", str(slideSize)], "WINDOW_PI")
    snpP, snpO, snpE = call_vcftools(inputVcf, ["--SNPdensity", str(windowSize)], "SNP_DENSITY")
    # missP, missO, missE = call_vcftools(inputVcf, ["--missing-site"], "MISSINGNESS")
    # frqP, frqO, frqE = call_vcftools(inputVcf, ["--freq"], "ALLELE_FREQUENCY")
    fullProcessList = [fstP, piP, snpP]
    fullProcessLogList = [fstO, fstE, piO, piE, snpO, snpE]
    print("\nDone initiating full VCF processes.")
    print("Checking if they are finished running...")
    for p in fullProcessList:
        if p.poll() == None:
            print("   Waiting on process...")
            p.wait()
    print("All full VCF processes are complete.")
    for logFile in fullProcessLogList:
        logFile.close()


    ### RUN STATS ON SEX-SPLIT VCFS
    
    #first check that splitting is complete
    check_wait(splitMaleProc, "splitMales")
    check_wait(splitFemaleProc, "splitFemales")
    print("\nSplitting VCF by sex is complete.")

    #START SPLIT BY SEX KEEPING ONLY VARIANT SITES FOR SNP DENSITY
    splitMaleSnpProc, maleSnpVcf, splitFemaleSnpProc, femaleSnpVcf, splitE, splitO = call_vcftools_split_sex_for_snp_density()


    #RUN VCF TOOLS ON SEX-SPECIFIC VERSIONS
    #start frequencies and site depth first
    #frequency
    ffrqP, ffrqO, ffrqE = call_vcftools(femaleVcf, ["--freq"], "ALLELE_FREQUENCY")
    mfrqP, mfrqO, mfrqE = call_vcftools(maleVcf, ["--freq"], "ALLELE_FREQUENCY")
    frqLogList = [ffrqO, ffrqE, mfrqO, mfrqE]
    #depth
    fdpthP, fdpthO, fdpthE = call_vcftools(femaleVcf, ["--site-depth"], "SITE_DEPTH")
    mdpthP, mdpthO, mdpthE = call_vcftools(maleVcf, ["--site-depth"], "SITE_DEPTH")



    #start other stats
    fpiP, fpiO, fpiE = call_vcftools(femaleVcf, ["--window-pi", str(windowSize), "--window-pi-step", str(slideSize)], "WINDOW_PI")
    mpiP, mpiO, mpiE = call_vcftools(maleVcf, ["--window-pi", str(windowSize), "--window-pi-step", str(slideSize)], "WINDOW_PI")
    fmissP, fmissO, fmissE = call_vcftools(femaleVcf, ["--missing-site"], "MISSINGNESS")
    mmissP, mmissO, mmissE = call_vcftools(maleVcf, ["--missing-site"], "MISSINGNESS")

    fdpthP, fdpthO, fdpthE = call_vcftools(femaleVcf, ["--site-depth"], "SITE_DEPTH")
    mdpthP, mdpthO, mdpthE = call_vcftools(maleVcf, ["--site-depth"], "SITE_DEPTH")





    #CHECK THAT THE SNP SPLITS ARE DONE, THEN GET SNP DENSITY FOR EACH SEX
    check_wait(splitMaleSnpProc, "splitMalesForSnps")
    check_wait(splitFemaleSnpProc, "splitFemalesForSnps")
    print("\nSplitting VCF by sex for SNPs is complete.")
    fsnpP, fsnpO, fsnpE = call_vcftools(femaleSnpVcf, ["--SNPdensity", str(windowSize)], "SNP_DENSITY")
    msnpP, msnpO, msnpE = call_vcftools(maleSnpVcf, ["--SNPdensity", str(windowSize)], "SNP_DENSITY")


    #get allele frequency statistics for the windows
    windowList = window_vcf(inputVcf, windowSize, slideSize, chrom, chromLength)

    #check that allele sex-specific allele frequency processes are complete
    check_wait(ffrqP, "female-allele-frequency")
    check_wait(mfrqP, "male-allele-frequency")
    check_wait(fdpthP, "female-depth")
    check_wait(mdpthP, "male-dpeth")
    print("\nGathering sex-specific allele frequencies is complete.")
    afRes = allele_frequency_stats(femaleVcf, maleVcf, windowList)


    ##break vcfs into windows
    p = Pool(Ncore)
    print("\nRunning window splitting with multiprocessing pool = {}:".format(Ncore))
    p.map(call_vcftools_window, range(len(windowList)))
    print("Done.")

    #ASSEMBLE THE OUTPUTS
    #check everything's done
    print('\nFinal check that all VCF tools processes are complete.')
    allProcList = [ffrqP, mfrqP, fpiP, mpiP, fmissP, mmissP, fsnpP, msnpP, fdpthP, mdpthP]
    for proc in allProcList:
        check_wait(proc, "...")
    #assemble outputs
    assemble_results()

    #REPORT TIME TO RUN
    print("\n\nFinished Running Functions.")
    print("chromosome = {}".format(chrom))
    print("window size = {}".format(windowSize))
    print("slide size = {}".format(slideSize))
    print("total windows = {}".format(len(windowList)))
    print("Assembled results saved as:  {}".format(finalOutput))
    print("Run Completed.")
    Time = time.time() - Start_time
    print('\nTime took to run: {}'.format(Time))



