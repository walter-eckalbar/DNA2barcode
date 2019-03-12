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
from functools import partial
path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path+'/modules')
from Bio import SeqIO
from barcodeExtract import *



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
barcodeRead = args.barcodeRead

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
print("barcode read is : " + bcRead)

# set up variables
alignmentBase = "alignOut"
alignments = alignmentBase+"/Aligned.sortedByCoord.out.bam"

#loading fasta file
faDict = SeqIO.to_dict(SeqIO.parse(library, "fasta"))
seqList = faDict.keys()

# create file in a directory of all reads aligning to a test sequence 
# then extract fastq reads from barcode file from those sequences

pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
func = partial(runExtractFastqByNames, bcRead, alignments, errorCorrectDir)
pool.map(func, seqList)
pool.close()
pool.join()	













