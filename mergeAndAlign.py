# mergeAndAlign.py #

# This script will merge then align the fastq files
# for assignment of BCs from pre-sequencing.

# input is directory of fastq files from pre-seq,
# the instert size of the libary (200bp) and the
# length of the sequencing reads (146bp). Assumes
# barcode read is R2


import subprocess, os, string
import csv
import math
import getopt
import sys
import os.path
import argparse
import re
path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+'/modules')
from flashMergeFastq import *
from starAlign import *

VERSION = 0.1

warningContinue = """WARNING: No arguments were specified, running with the following defaults:
library=library/reference.fa
threads=4
fastqDir=fastq
insertSize=200
readLength=146
barcodeRead=R2
Are you sure you want to continue? (Y/N): """

library = ''
threads = ''
fastqDir = ''
insertSize = ''
readLength = ''
barcodeRead= ''

parser = argparse.ArgumentParser(description="Align fastq reads.")
parser.add_argument("-l", "--library", default="library/reference.fa", type=str, help="-l [--library], Choose fasta file to great reference from")
parser.add_argument("-f", "--fastqDir", default="fastq", type=str, help="-f [--fastq_dir], Choose fastq directory for merging and alignment")
parser.add_argument("-i", "--insertSize",default="200",type=int, help="-i [--insertSize], Length of test sequences. Default 200bp")
parser.add_argument("-r", "--readLength", default="146",type=int, help="-r, [--readLength], Length of paired end reads. Default 146bp")
parser.add_argument("-t", "--threads", default=4, type=int, help="-t [--threads], genomeGenerate on this many threads.")
parser.add_argument("-bc", "--barcodeRead", default="R2",type=str, help="-bc [--barcodeRead], read with enhancer barcode")
parser.add_argument("--version", action="store_true", help="[--version], prints the version of the script and then exits")
args = parser.parse_args()

if args.version == True:
	print "version: "+ str(VERSION)
	print "exiting..."
	sys.exit()

library = args.library
threads = args.threads
fastqDir = args.fastqDir
insertSize = args.insertSize
readLength = args.readLength
barcodeRead = args.barcodeRead

print("selected library file is : " + library)
print("number of threads is : " + str(threads))

# Setting up libraries #
STAR_dir = library.split("/")[1]
STAR_dir = STAR_dir.split(".fa")[0]
STAR_dir = "STAR_" + STAR_dir


# Some math on the max and min overlap
overlap = 2*readLength - insertSize
print("theoretical overlap is : " + str(overlap))
maxOverlap = overlap + 2
minOverlap = overlap - 2 # doing this because the reads seem to slip a bit and we often get n - 1 overlap, n + 1 seems rare

# Running overlap
try:
	targets = os.listdir(fastqDir)
except OSError:
	print("no fastq directory, exiting for now.... (this should be fixed for countings step)")
	print "cur fastqDir: " + fastqDir
	sys.exit()
for i in targets:
	if i.find("_R1") > 1:
		R1 = i
	if i.find("_"+barcodeRead) > 1:
		targets.remove(i)
	continue

fastqMergeBase = fastqDir+"_merge"
try: 
	os.mkdir(fastqMergeBase)
except OSError:
	pass

mergeOutputBase =  fastqMergeBase+"/"+R1.split("_R1_")[0]

print("Merged file basename will be : " + mergeOutputBase)
print("Merging files : " + targets[0] + " and " + targets[1])

targetsFinal = [fastqDir + "/" + x for x in targets]

runOverlap(targetsFinal,minOverlap,maxOverlap,mergeOutputBase,threads)


# Run STAR alignment from merged fastqs
alignmentBase = "alignOut"
try: 
	os.mkdir(alignmentBase)
except OSError:
	pass

mergedFastq = mergeOutputBase+".extendedFrags.fastq.gz"

runSTAR(mergedFastq, STAR_dir, threads, alignmentBase)




