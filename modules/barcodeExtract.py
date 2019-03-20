#!/usr/bin/env python

# This script will overlap the PE reads into single merged read

import sys, re, math, subprocess
import os
import pysam
from Bio import Align
from Bio import SeqIO

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


def runFilterWindowBCs(errorCorrectDir, minBarcodeCounts,percentBCReads, window):
	j = window
	tmpDict = {}
	aboveMinCutoffDict = {}
	indFastq = "readsOn."+j+".fastq"
	ecFastq = errorCorrectDir+"/"+j+"/"+"readsOn."+j+".ec.fastq"
	if os.path.isfile(ecFastq):
		pass
	else:
		#print "ln -s " + indFastq + " " + ecFastq
		subprocess.call("ln -s " + indFastq + " " + ecFastq, shell=True)

	fastqDict = SeqIO.to_dict(SeqIO.parse(ecFastq,"fastq"))
	
	counter = 0
	for key in fastqDict:
		counter += 1 
		seqK = fastqDict[key].seq
		if seqK in tmpDict:
			tmpDict[seqK] = tmpDict[seqK] + 1
		else:
			tmpDict[seqK] = 1

	sumReads = len(fastqDict)
	sumReadFilteredReads = sumReads
	
	if bool(tmpDict):
		valueMax = max(tmpDict.itervalues())
	
		iCoverage = 1
		for iCoverage in range(1,valueMax):
			delList = []
			tmp2Dict = dict( (k, v) for k, v in tmpDict.items() if v >= iCoverage)
			sumReadFilteredReads = sum(tmp2Dict.itervalues())
			if (sumReadFilteredReads <= (sumReads * percentBCReads)):
				break
			else:
				continue
		#try:
		#	tmp2Dict
		#	if isinstance(tmp2Dict,dict):# in locals(): #if isinstance(ele,dict)
		for tmp2Key in tmp2Dict:
			if tmp2Dict[tmp2Key] >= minBarcodeCounts:
				aboveMinCutoffDict[tmp2Key] = tmp2Dict[tmp2Key]
			
			for tmp3Key in aboveMinCutoffDict:
				tmp3List = [j, aboveMinCutoffDict[tmp3Key]]
				return tmp3Key, tmp3List
					#return tmp3Key, tmp3List
					#if tmp3Key in finalDict:
					#	finalDict[tmp3Key].append(tmp3List)
					#else:
					#	finalDict[tmp3Key] = list()
					#	finalDict[tmp3Key].append(tmp3List)
#		except UnboundLocalError:
#			pass








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