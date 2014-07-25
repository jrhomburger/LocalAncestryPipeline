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

splitchrs = parse_param("splitchrs")
if splitchrs == "False"
	splitchrs = False
else:
	splitchrs = True

remake_map = parse_param("remake_map")
if remake_map == "False":
	remake_map = False
else:
	remake_map = True

rfmix_run = parse_param("rfmix_run")
if rfmix_run == "False":
	rfmix_run = False
else:
	rfmix_run = True

logfolder = parse_param("logfolder", allparams)

if shapeitphase:
	phased_adm = shapeItPhase(plinkbfile, genetic_map, "", max_threads=max_threads)
	phased_refs = []
	for k in ref_plinkfiles:
		phased_refs.append(shapeItPhase(k + "_chr" + chr, genetic_map, "", max_threads=max_threads))

if rmfixrun:
	rfmix_in = shapeItToRFMixMultiClass(phased_adm, phased_refs)

if remake_map:
	this_map = makeSNPMap(plinkbfile + ".bim", genetic_map)

rfmix_out = runRFMixRephase(rfmix_in[0], rfmix_in[1], this_map, window_size, generations, em_iters, max_threads=max_threads)

# Write tempfile to communicate with main script 
tempfile = open(logfolder + "/" + "out_chr" + chr, "w")
tempfile.write(rfmix_out)
tempfile.write("\n")
tempfile.close()

