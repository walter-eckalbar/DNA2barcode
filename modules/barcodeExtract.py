#!/usr/bin/env python

# This script will overlap the PE reads into single merged read

import sys, re, math, subprocess
import os
import pysam
from Bio import Align

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
	ecOutputErr = outputBase+".ec.err.txt"
	ecOutputLog = outputBase+".ec.log.txt"

	#print("seqkit grep --pattern-file "+ readListFile + " " + fastq + " > " + outputFastq + " 2> " + outputBase+".err.txt")
	subprocess.call("seqkit grep --pattern-file "+ readListFile + " " + fastq + " > " + outputFastq + " 2> " + outputBase+".err.txt", shell=True)
	#subprocess.call("ace 1500 " + outputFastq + " " + ecOutput + " > " + ecOutputLog + " 2> " + ecOutputErr,shell=True)



def alignSeqs(keySequences, editDistance, count):
	aligner = Align.PairwiseAligner()
	aligner.open_gap_score=-0.5
	aligner.extend_gap_score=-0.5
	minScore = 15 - editDistance
	i = keySequences[int(count)]
	countj = int(int(count)+1)
	for j in keySequences[countj:]:
		alignments = aligner.align(i,j)
		if alignments.score >= minScore and  alignments.score < 15:
			alignment = alignments[0]
			return i,j