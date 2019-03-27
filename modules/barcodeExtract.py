#!/usr/bin/env python

# This script will overlap the PE reads into single merged read

import sys, re, math, subprocess
import os,csv
import pysam
from Bio import Align
from Bio import SeqIO

def findSumThreshold(vect,perc):
    v = sorted(vect, reverse=True)
    sumV = sum(v)
    vFiltered =[]
    for i in v:
        vFiltered = [j for j in v if j >= i]
        sumVFiltered = sum(vFiltered)
        if sumVFiltered >= (sumV * float(perc)):
            return i
        else:
            continue

def getUniqBCsPerWindow(errorCorrectDir,window):
	j = window
	tmpDict = {}
	
	uniqueFile = errorCorrectDir+"/"+j+"/"+"uniqBarcodes.counts."+j+".txt"
	indFastq = "readsOn."+j+".fastq"
	ecFastq = errorCorrectDir+"/"+j+"/"+"readsOn."+j+".ec.fastq"
	if os.path.isfile(ecFastq):
		pass
	else:
		subprocess.call("ln -s " + indFastq + " " + ecFastq, shell=True)

	fastqDict = SeqIO.to_dict(SeqIO.parse(ecFastq,"fastq"))
	
	counter = 0
	for key in fastqDict:
		counter += 1 
		seqK = fastqDict[key].seq
		if "N" not in seqK:
			if seqK in tmpDict:
				tmpDict[seqK] = tmpDict[seqK] + 1
			else:
				tmpDict[seqK] = 1
		else:
			continue

	uniqueFileOpen=open(uniqueFile,'wb')
	uniqueFileExp=csv.writer(uniqueFileOpen, lineterminator='\n', delimiter = '\t')
	uniqueFileExp.writerow(["Window", "Barcode Seq", "counts"])
	for key in tmpDict:
		value = tmpDict[key]
		uniqueFileExp.writerow([j, key, value])
	return tmpDict


def allUniqBCbyWindow(errorCorrectDir,windowList):
	uniqBCdict = {}
	for i in windowList:
		uniqueFile = errorCorrectDir+"/"+i+"/"+"uniqBarcodes.counts."+i+".txt"
		with open(uniqueFile) as openFile:
			for line in openFile.readlines():
				line = line.split('\n')[0]
				line = line.split('\t')

				
				window = line[0]
				barcode = line[1]
				coverage = line[2]
				key = window + "_" + barcode
				if window != "Window":
					uniqBCdict[key] = [window,barcode,coverage]
    	#uniqueFile.close()
	
	uniqueAllFile = "final_outputs/uniqBarcodesAll.counts.txt"
	uniqueAllFileOpen=open(uniqueAllFile,'wb')
	uniqueAllFileExp=csv.writer(uniqueAllFileOpen, lineterminator='\n', delimiter = '\t')
	uniqueAllFileExp.writerow(["Window_BC","Window", "Barcode Seq", "counts"])
	for key in uniqBCdict:
		value = uniqBCdict[key]
		window = value[0]
		barcode = value[1]
		counts = value[2]
		uniqueAllFileExp.writerow([key,window,barcode,counts])




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
	#subprocess.call("seqkit grep --pattern-file "+ readListFile + " " + fastq + " > " + outputFastq + " 2> " + outputBase+".err.txt", shell=True)
	subprocess.call("samtools fqidx -r "+ readListFile + " -o " + outputFastq + " " + fastq + " > "  + outputBase+".log.txt" + " 2> " + outputBase+".err.txt", shell=True)
	subprocess.call("ace 1500 " + outputFastq + " " + ecOutput + " > " + ecOutputLog + " 2> " + ecOutputErr,shell=True)


def runFilterWindowBCs(errorCorrectDir, minBarcodeCounts,percentBCReads, window):
	j = window
	tmpDict = {}
	aboveMinCutoffDict = {}
	indFastq = "readsOn."+j+".fastq"
	ecFastq = errorCorrectDir+"/"+j+"/"+"readsOn."+j+".ec.fastq"

	tmpDict = getUniqBCsPerWindow(errorCorrectDir, window)

	
	if bool(tmpDict):
		coverageList = list(tmpDict.values())
		percLimit = findSumThreshold(coverageList, percentBCReads)


		tmp2Dict = {}

		for k in tmpDict:
			v = tmpDict[k]
			if v >= percLimit:
				tmp2Dict[k] = v

		outputList = []
		tmp3Dict = {}
		for tmp2Key in tmp2Dict:
			if tmp2Dict[tmp2Key] >= minBarcodeCounts:
				tmp2List = [j, tmp2Dict[tmp2Key]]
				tmp3Dict[tmp2Key] = tmp2List
		for tmp3Key in tmp3Dict:
			tmp3Value = tmp3Dict[tmp3Key]
			outputList.append([tmp3Key,tmp3Value])
			#return tmp3Key, tmp3Value
		return outputList





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