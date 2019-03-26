#!/usr/bin/env python

# This script will overlap the PE reads into single merged read

import sys, re, math, subprocess, os

def runOverlap(fastqList, minOverlap, maxOverlap, outPutName, threads):
	R1 = fastqList[0]
	R2 =  fastqList[1]
	if os.path.exists(outPutName + ".err") and os.path.getsize(outPutName + ".err") == 0:
		print("skipping fastq merge, file exists")
	else:
		print("flash  -o  " + outPutName + " -z -t " + str(threads) + " -M " + str(maxOverlap) + " -m " + str(minOverlap) + " -x .33 " + R1 + " " + R2)
		subprocess.call("flash  -o  " + outPutName + " -z -t " + str(threads) + " -M " + str(maxOverlap) + 
			" -m " + str(minOverlap) + " -x .33 " + R1 + " " + R2 + " > " + outPutName + ".log 2> " + outPutName + ".err", shell=True)

		