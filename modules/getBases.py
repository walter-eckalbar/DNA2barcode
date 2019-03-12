#!/usr/bin/env python

# This script will identify some genomeGenerate STAR parameters

import sys, re, math


def getSTARGenomeParameters(inFile):
	GENOMELENGTH = 3000000000
	DNApattern = re.compile("^[ACGTU]*$")
	reads = 0
	allChars = ""

	StarParams = []
	with open(inFile, 'rU') as infile:
	    for line in infile:
	    	if DNApattern.search(line.rstrip("\n")):
	    		allChars+=line.rstrip("\n")
	    	if '>' in line.rstrip("\n"):
	    		reads+=1
	
	print "num bases: "+"{:,}".format(len(allChars))
	print "num reads: "+"{:,}".format(reads)
	genomeSAindexNbasesN = min(14, math.log(len(allChars), 2)/2-1)
	genomeChrBinNbitsN = min(18, math.log(len(allChars)/reads, 2))
	StarParams.append(genomeSAindexNbasesN)
	StarParams.append(genomeChrBinNbitsN)
	
	print "suggested genomeSAindexNbases value: "+str( min(14, math.log(len(allChars), 2)/2-1) ) #min(14, log2(GenomeLength)/2 - 1)
	print "suggested genomeChrBinNbits value: "+str( min(18, math.log(len(allChars)/reads, 2)) ) #min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])
	return StarParams