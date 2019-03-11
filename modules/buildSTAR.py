#!/usr/bin/env python

# This script will build the STAR index files

import sys, re, math, subprocess

def makeStarIndex(library, threads, outDir, StarParams):
	SAindexNbases = int(round(StarParams[0]))
	ChrBinNbits =  int(round(StarParams[1]))
	print("STAR --runMode genomeGenerate --runThreadN " + str(threads) + 
		" --genomeDir " + outDir + " --genomeFastaFiles " + library + 
		" --genomeSAindexNbases " + str(SAindexNbases) + " --genomeChrBinNbits " + str(ChrBinNbits))
	subprocess.call("STAR --runMode genomeGenerate --runThreadN " + str(threads) + 
		" --genomeDir " + outDir + " --genomeFastaFiles " + library + 
		" --genomeSAindexNbases " + str(SAindexNbases) + " --genomeChrBinNbits " + str(ChrBinNbits), shell=True)