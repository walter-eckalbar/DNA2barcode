#!/usr/bin/env python

# This script will overlap the PE reads into single merged read

import sys, re, math, subprocess
import os
import pysam

def runExtractFastqByNames(fastq, bam, errorCorrectDir, seqID):

	readListFile = errorCorrectDir+"/"+seqID+"/"+"readsOn."+seqID+".txt"
	outputFastq = errorCorrectDir+"/"+seqID+"/"+"readsOn."+seqID+".fastq"
	
	#loading alignments
	bamFile = pysam.AlignmentFile(bam, "rb")

	tmpList = []
	try: 
		os.mkdir(errorCorrectDir+"/"+seqID)
	except OSError:
		pass
	tmpSam = bamFile.fetch(seqID)
	for x in tmpSam:
		tmpList.append(x.query_name)
	with open(readListFile, 'w') as f:
		for item in tmpList:
			f.write("%s\n" % item)

	outputBase = outputFastq.split(".fastq")[0]
	ecOutput = outputBase+".ec.fastq"
	ecOutpubtErr = outputBase+".ec.err.txt"
	ecOutpubtLog = outputBase+".ec.log.txt"

	#print("seqkit grep --pattern-file "+ readListFile + " " + fastq + " > " + outputFastq + " 2> " + outputBase+".err.txt")
	subprocess.call("seqkit grep --pattern-file "+ readListFile + " " + fastq + " > " + outputFastq + " 2> " + outputBase+".err.txt", shell=True)
	subprocess.call("ace 1500 " + outputFastq + " " + ecOutput + " > " + ecOutpubtLog + " 2> " + ecOutputErr,shell=True)