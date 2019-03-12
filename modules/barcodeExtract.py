#!/usr/bin/env python

# This script will overlap the PE reads into single merged read

import sys, re, math, subprocess

def runExtractFastqByNames(fastq, nameList, output):
	#print("seqkit grep --pattern-file "+ nameList + " " + fastq + " > " + output)
	outputBase = output.split(".fastq.gz")[0]
	subprocess.call("seqkit grep --pattern-file "+ nameList + " " + fastq + " | gzip -c > " + output " 2> " + outputBase+".err.txt", shell=True)