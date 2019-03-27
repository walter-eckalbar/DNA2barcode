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
import pysam
import multiprocessing
from multiprocessing import Manager,Process
import gzip
import time
import datetime
import gc
from functools import partial
path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+'/modules')
from Bio import SeqIO
from Bio import Align
import itertools
from collections import defaultdict
from barcodeExtract import *



VERSION = 0.1

warningContinue = """WARNING: No arguments were specified, running with the following defaults:
library=library/reference.fa
threads=4
fastqDir=fastq
insertSize=200
readLength=146
barcodeRead=R2
minBarcodeCounts=20
editDistanceMax=2
Are you sure you want to continue? (Y/N): """

def getTime():
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	return(st)

library = ''
threads = ''
fastqDir = ''
insertSize = ''
readLength = ''
minBarcodeCounts = ''
barcodeRead= ''
editDistanceMax= ''
percentBCReads= ''

parser = argparse.ArgumentParser(description="Align fastq reads.")
parser.add_argument("-l", "--library", default="library/reference.fa", type=str, help="-l [--library], Choose fasta file to great reference from")
parser.add_argument("-f", "--fastqDir", default="fastq", type=str, help="-f [--fastq_dir], Choose fastq directory for merging and alignment")
parser.add_argument("-t", "--threads", default=4, type=int, help="-t [--threads], genomeGenerate on this many threads (default : 4).")
parser.add_argument("-bc", "--barcodeRead", default="R2",type=str, help="-bc [--barcodeRead], read with enhancer barcode (default : R2)")
parser.add_argument("-m", "--minBarcodeCounts", default="20",type=int, help="-m [--minBarcodeCounts], minimum reads of barcode to keep\n (default : 20)")
parser.add_argument("-e", "--editDistanceMax", default="2",type=int, help="-e [--editDistanceMax], maximum edit distance between white list barcodes\n default  : 2")
parser.add_argument("-p", "--percentBCReads", default="92.5",type=float, help="-p [--percentBCReads], BCs that contain the top X (default 92.5) percentile of reads \n(removes many frequent low count BCs)")
parser.add_argument("--version", action="store_true", help="[--version], prints the version of the script and then exits")
args = parser.parse_args()

if args.version == True:
	print "version: "+ str(VERSION)
	print "exiting..."
	sys.exit()

library = args.library
threads = args.threads
fastqDir = args.fastqDir
fastqDir = fastqDir.split("/")[0]
barcodeRead = args.barcodeRead
minBarcodeCounts = args.minBarcodeCounts
editDistanceMax = args.editDistanceMax
percentBCReads = args.percentBCReads
percentBCReads = percentBCReads/100

print("selected library file is : " + library)
print("number of threads is : " + str(threads))

# Setting up libraries #
STAR_dir = library.split("/")[1]
STAR_dir = STAR_dir.split(".fa")[0]
STAR_dir = "STAR_" + STAR_dir

errorCorrectDir = "errorCorrect"
try: 
	os.mkdir(errorCorrectDir)
except OSError:
	pass

finalDir = "final_outputs"
try: 
	os.mkdir(finalDir)
except OSError:
	pass

#identify bc read
try:
	targets = os.listdir(fastqDir)
except OSError:
	print("no fastq directory, exiting for now.... (this should be fixed for countings step)")
	print "cur fastqDir: " + fastqDir
	sys.exit()
for i in targets:
	if i.find("_"+barcodeRead) > 1:
		bcRead = i
	continue
bcRead = fastqDir+"/"+bcRead
bcReadFastq = bcRead.split(".gz")[0]
print("barcode read is : " + bcRead)


if os.path.isfile(bcReadFastq):
	print("unzipped file exists, skipping")
else:
	print("gunzip-ing fastq : " + getTime())
	subprocess.call("unpigz -c -p " + str(threads) + " " + bcRead + " > " + bcReadFastq,shell=True)

bcReadFastqIDX = bcReadFastq + ".fai" 

if os.path.isfile(bcReadFastqIDX):
	print("fastq index file exists, skipping")
else:
	print("indexing fastq : " + getTime())
	subprocess.call("samtools fqidx " + bcRead + " > " + bcReadFastqIDX+".log.txt 2> " + bcReadFastqIDX+".err.txt",shell=True)

# set up variables
alignmentBase = "alignOut"
alignments = alignmentBase+"/Aligned.sortedByCoord.out.bam"

#loading fasta file
faDict = SeqIO.to_dict(SeqIO.parse(library, "fasta"))
windowList = faDict.keys()

del faDict
gc.collect()


# create file in a directory of all reads aligning to a test sequence 
# then extract fastq reads from barcode file from those sequences
print("running multiprocessing of barcode error correction : " + getTime())
pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
func = partial(runExtractFastqByNames, bcRead, alignments, errorCorrectDir)
pool.map(func, windowList)
pool.close()
pool.join()

# read error corrected reads, create dictionary of library seqs, values of BCs and counts
# dictionary of keys for Windows, values are a second dictionary with seqs as keys and counts as values 
print("getting BC key, window value dictionaries : " + getTime())
bcDictSeqsAsKeys = {}

def runMultiThreadedReading():
	pool = multiprocessing.Pool(int(threads))
	manager = Manager()
	#bcDictSeqsAsKeys = manager.dict()
	#func = partial(runFilterWindowBCs, bcDictSeqsAsKeys,errorCorrectDir, minBarcodeCounts, percentBCReads)
	func = partial(runFilterWindowBCs, errorCorrectDir, minBarcodeCounts, percentBCReads)
	listSeqs = pool.map(func, windowList)

	listSeqsNoneRM = []
	for i in listSeqs:
		if i != None:
			listSeqsNoneRM.append(i)
	
	flat_listSeqs = [item for sublist in listSeqsNoneRM for item in sublist]


	for i in flat_listSeqs:
		if i != None:
			key = i[0]
			value = i[1]
			if key in bcDictSeqsAsKeys:
				bcDictSeqsAsKeys[key].append(value)
			else:
				bcDictSeqsAsKeys[key] = [value]

	allUniqBCbyWindow(errorCorrectDir,windowList)


runMultiThreadedReading()


blackList = []
for seqK in bcDictSeqsAsKeys:
	if len(bcDictSeqsAsKeys[seqK]) > 1:
		blackList.append(seqK)

print("number of blackListed barcodes : " + str(len(blackList)))

bcDictSeqsAsKeys
for seqK in blackList:
	del bcDictSeqsAsKeys[seqK]

print("number of barcodes post 1st round blacklist filtering : " + str(len(bcDictSeqsAsKeys)))


# loop through bcs and remove those with small edit distances to a second read.
editDistBlackList = []

print("aligning BCs to eachother :" + getTime())

keySeqs = bcDictSeqsAsKeys.keys()
seqRanges = range(0,len(keySeqs))


pool = multiprocessing.Pool(int(threads))
func = partial(alignSeqs, keySeqs,editDistanceMax)
alignedBlackList = pool.map(func, seqRanges)
pool.close()
pool.join()

print("Final barcode whitelisting :" + getTime())

alignedBlackList  = [x for x in alignedBlackList if x is not None]

# Flatten the list of lists
alignedBlackListFlat = [item for sublist in alignedBlackList for item in sublist]
alignedBlackListUniq = list(set(alignedBlackListFlat))

print("number of barcodes too close between windows : " + str(len(alignedBlackListUniq)))

#removing black listed bc based on edit distance
for iEditBlackList in alignedBlackListUniq:
	del bcDictSeqsAsKeys[iEditBlackList]

print("number of barcodes post 2nd round blacklist filtering : " + str(len(bcDictSeqsAsKeys)))

fastaFileOut=open(finalDir+"/finalBarcodes.fa",'wb')
tableFileOut = open(finalDir+"/finalBarcodes.counts.txt", 'wb')

tableFileExp=csv.writer(tableFileOut, lineterminator='\n', delimiter = '\t')
fastaFileExp=csv.writer(fastaFileOut, lineterminator='\n', delimiter = '\t')

tableFileExp.writerow(["test_sequences","barcode_ID","barcode_seq","barcode_count"])

for key in bcDictSeqsAsKeys:
	bcSeq = key
	test_seq = bcDictSeqsAsKeys[key]
	windowAndCoverage = test_seq[0]
	window = windowAndCoverage[0]
	bcCoverage = windowAndCoverage[1]
	barcode_id = window+"_"+bcSeq
	tableFileExp.writerow([window , barcode_id,bcSeq,str(bcCoverage)])
	fastaFileExp.writerow([">"+barcode_id])
	fastaFileExp.writerow([bcSeq])


print("All finished :" + getTime())




















