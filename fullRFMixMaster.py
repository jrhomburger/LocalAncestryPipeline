### This is the master program for the RFMix pipeline.
### Its job is to call the other RFMix pipeline programs
### This could be either through qsub or not through qsub

import sys
import re
import os
import os.path
import time
from rfmixPipelineFunctions import *

paramfile = open(sys.argv[1], "r")
allparams = paramfile.read()

## Param qsub, default is False. If True, will qsub
qsub = parse_param("qsub", allparams)
if (qsub == "True"):
	qsub = True
else:
	qsub = False

max_threads = parse_param("maxthreads", allparams)
if max_threads == "":
	max_threads = 1

phasemix = parse_param("phasemix", allparams)
if (phasemix == "False"):
	phasemix = False
else:
	phasemix = True


run_concat = parse_param("runconcat", allparams)
if (run_concat == "False"):
	run_concat = False
else:
	run_concat = True

requality = parse_param("requality", allparams)
if requality == "True":
	run_concat=False
	phasemix = False
	requality = True
else:
	requality = False

## First take as input the parameter file

splitchrs = parse_param("splitchrs")
if splitchrs == "False"
	splitchrs = False
else:
	splitchrs = True


run_phase_mix=phasemix

chroms = parse_param("chroms", allparams)
if chroms == "":
	chroms = range(1,23)
else:
	chroms = chroms.split(",")
print chroms
# Test on these two chromosomes locally

logfolder = parse_param("logfolder", allparams)
if not os.path.exists(logfolder) and logfolder != "":
	os.makedirs(logfolder)
	logfolder = logfolder + "/"


## Split into chromosomes:
print allparams
plinkbfile = parse_param("plinkbfile", allparams)
#print plinkbfile
ref_plinkfiles = parse_ref_files("mixreference", allparams)
#print ref_plinkfiles

all_plink = ref_plinkfiles + [plinkbfile]
#print all_plink
if run_phase_mix:
	if splitchrs: # Split plink files into chromosomes
		for thisfile in all_plink:
			for chrom in chroms:
				out = thisfile + "_chr" + str(chrom)
				os.system("plink --noweb --make-bed --bfile " + thisfile + " --chr " + str(chrom) + " --out " + out)

	## Run each chromosome file:
	#print "there"
	for chrom in chroms:
		if qsub:
			logfile = logfolder + "mix" + str(chrom)
			name = "mix" + str(chrom)
			print("qsub -V -cwd -v par=" + sys.argv[1] + ",chr=" + str(chrom) + " -o " + logfile + " -e " + logfile +
				" -N " + name + " -pe shm " + max_threads + " singleChrRFMix.sh")
			os.system("qsub -V -cwd -v par=" + sys.argv[1] + ",chr=" + str(chrom) + " -o " + logfile + " -e " + logfile +
				" -N " + name + " -q extended -pe shm " + max_threads + " singleChrRFMix.sh")
		else:
			os.system("./singleChrRFMix.sh " + sys.argv[1] + " " + str(chrom))

	running = [True] * len(chroms)
	while any(running):
		# Check for existence of all output files:
		for k, chrom in enumerate(chroms):
			if (os.path.exists(logfolder + "out_chr" + str(chrom))): # check which ones aren't running
				running[k] = False
		if any(running):
			time.sleep(60) # wait a minute if some are running

	# delete temp files at this point 
	for chrom in chroms:
		os.remove(logfolder + "out_chr" + str(chrom))


## Concatenate output and format for PCAmask:

if run_concat:
	if qsub:
		print "qsub"
		logfile = logfolder + "/concat"
		name = "concat"
		os.system("qsub -V -cwd -v par=" + sys.argv[1] + " -o " + logfile + " -e " + logfile +
				" -N " + name + " -l h_vmem=24G run_concat.sh")
	else:
		os.system("./run_concat.sh " + sys.argv[1])

	run = True
	while run: # wait for concat to finish
		if os.path.exists( plinkbfile + "_cat"):
			run = False
		if run:
			time.sleep(60)
	os.remove(plinkbfile + "_cat")		
		
if requality:
	if qsub:
		logfile = logfolder + "requality"
		os.system("qsub -V -cwd -v par=" + sys.argv[1] + " -o " + logfile + " -e " + logfile +
				" -N requality -l h_vmem=10G" + " requality.sh")
	else:
		os.system("./requality.sh " + sys.argv[1])
