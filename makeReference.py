# makeReference.py #

# This script will create STAR reference for the 
# assignment of BCs from pre-sequencing.

# input is fasta file of library sequences 

path = "/Users/waltereckalbar/tools/MPRAssign/"


import subprocess, os, string
import csv
import math
import getopt
import sys
import os.path
import argparse
import re
sys.path.insert(0, path+'modules')
from getBases import *
from buildSTAR import *

VERSION = 0.1

warningContinue = """WARNING: No arguments were specified, running with the following defaults:
library=library/reference.fa
threads=4
Are you sure you want to continue? (Y/N): """

library = ''

parser = argparse.ArgumentParser(description="Align fastq reads.")
parser.add_argument("-l", "--library", default="library/reference.fa", type=str, help="-l [--library], Chose fasta file to great reference from")
parser.add_argument("-t", "--threads", default=6, type=int, help="-t [--threads], genomeGenerate on this many threads.")
parser.add_argument("--version", action="store_true", help="[--version], prints the version of the script and then exits")
args = parser.parse_args()

if args.version == True:
	print "version: "+ str(VERSION)
	print "exiting..."
	sys.exit()

library = args.library
threads = args.threads

print("selected library file is : " + library)
print("number of threads is : " + str(threads))

# Determining some STAR parameters #
StarParams = getSTARGenomeParameters(library)

# Setting up libraries #
STAR_dir = library.split("/")[1]
STAR_dir = STAR_dir.split(".fa")[0]
STAR_dir = "STAR_" + STAR_dir

print("making STAR index directory : " + STAR_dir)
try: 
	os.mkdir(STAR_dir)
except OSError:
	pass

# Running STAR genomeGenerate command #
if os.path.isfile(STAR_dir+"/"+"SA"):
	print("STAR index already made, skipping")
else:
	makeStarIndex(library, threads, STAR_dir, StarParams)




