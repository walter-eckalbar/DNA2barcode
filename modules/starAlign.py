#!/usr/bin/env python

# This script will align the overlapped PE reads using STAR

import sys, re, math, subprocess

def runSTAR(fastq, starDir, threads, outPutName):
	parameters = " --alignIntronMin 20 --alignIntronMax 19 --alignMatesGapMax 10 --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate"
	parameters = parameters + " --outFilterType Normal --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNmax 999 --outFilterMultimapNmax 100"
	parameters = parameters + " --limitBAMsortRAM 23052597265 --outFilterMultimapScoreRange 1"
	print("STAR --runThreadN " + str(threads) + " --genomeDir " + starDir + " --readFilesIn " + fastq + " --outFileNamePrefix " + outPutName+"/" + parameters)
	subprocess.call("STAR --runThreadN " + str(threads) + " --genomeDir " + starDir + " --readFilesIn " + fastq + " --outFileNamePrefix " + outPutName+"/" + parameters,shell=True)
	subprocess.call("samtools index " + outPutName+"/Aligned.sortedByCoord.out.bam", shell=True)