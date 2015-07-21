### This script runs a single chromosome through RFMix
### It takes as arguments the original parameter file and the chromosome number

### By Julian Homburger

import sys
import re
from rfmixPipelineFunctions import * # import my functions

# Parameter file and chromosome on first line
paramfile = open(sys.argv[1], "r")
chr = sys.argv[2]
allparams = paramfile.read()

plinkbfile = parse_param("plinkbfile", allparams) + "_chr" + chr
genetic_map = re.sub(r"NNN", chr, parse_param("geneticmap", allparams)) # specify genetic map as type genetic_map_chrNNN_b36.txt
# will automatically replace NNN with necessary chromosome number

ref_plinkfiles = parse_ref_files("mixreference", allparams)

window_size = parse_param("windowsize", allparams)
if window_size == "":
	window_size = 0.2
	
generations = parse_param("generations", allparams)
if generations == "":
	generations = 8

em_iters = parse_param("emiterations", allparams)
if em_iters == "":
	em_iters = 2

max_threads = parse_param("maxthreads", allparams)
if max_threads == "":
	max_threads = 1

splitchrs = parse_param("splitchrs", allparams)
if splitchrs == "False":
	splitchrs = False
else:
	splitchrs = True


shapeitphase = parse_param("shapeitphase", allparams)
if shapeitphase == "False":
	shapeitphase = False
else:
	shapeitphase = True

remake_map = parse_param("remakemap", allparams)
if remake_map == "False":
	remake_map = False
else:
	remake_map = True

rfmix_run = parse_param("rfmixrun", allparams)
if rfmix_run == "False":
	rfmix_run = False
else:
	rfmix_run = True
print rfmix_run

duohmm = parse_param("duohmm", allparams)
if duohmm == "True":
	duohmm = True
else:
	duohmm = False
print duohmm

trioPhase = parse_param("triophase", allparams)
if trioPhase == "True":
	trioPhase = True
else:
	trioPhase = False
if (trioPhase):
	print "Using trio phase RFMix"
else:
	print "Using pop phase RFMix"


logfolder = parse_param("logfolder", allparams)

if shapeitphase:
	#print "here"
	phased_adm = shapeItPhase(plinkbfile, genetic_map, "", max_threads=max_threads, duohmm=duohmm)
	phased_refs = []
	for k in ref_plinkfiles:
		phased_refs.append(shapeItPhase(k + "_chr" + chr, genetic_map, "", max_threads=max_threads, duohmm=duohmm))
		
	rfmix_in = shapeItToRFMixMultiClass(phased_adm, phased_refs)
else:
	phased_adm = plinkbfile + "_shapeout"
	rfmix_in = [phased_adm + ".alleles", phased_adm + ".classes", phased_adm + ".amaps"]  


	
if remake_map:
	this_map = makeSNPMap(plinkbfile + ".bim", genetic_map)

if rfmix_run:
	rfmix_out = runRFMixRephase(rfmix_in[0], rfmix_in[1], this_map, window_size, generations, em_iters, max_threads=max_threads, trioPhase=trioPhase)

# Write tempfile to communicate with main script 
tempfile = open(logfolder + "/" + "out_chr" + chr, "w")
tempfile.write("Finished")
tempfile.write("\n")
tempfile.close()

