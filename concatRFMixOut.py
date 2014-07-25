### Program that concatenates RFMix Output and generates beagle input files 
### for input into PCAMask


import sys
import re
from rfmixPipelineFunctions import *

paramfile = open(sys.argv[1], "r")
allparams = paramfile.read()

chroms = parse_param("chroms", allparams)
if chroms == "":
	chroms = range(1,23)
else:
	chroms = chroms.split(",")

plinkbfile_pref = parse_param("plinkbfile", allparams)
window_size = parse_param("windowsize", allparams)
if window_size == "":
	window_size = 0.2
window_size = float(window_size)
	
generations = parse_param("generations", allparams)
if generations == "":
	generations = 8
generations = int(generations)

em_iters = parse_param("emiterations", allparams)
if em_iters == "":
	em_iters = 2
em_iters = int(em_iters)

qual_cut = parse_param("qualitycutoff", allparams)
if qual_cut == "":
	qual_cut = 0.95
qual_cut = float(qual_cut)

anc_cut = parse_param("ancestrycutoff", allparams)
if anc_cut == "":
	anc_cut = 0.25
else:
	anc_cut = float(anc_cut)

exclude_inds = parse_param("excludeinds", allparams)
exclude_sites = parse_param("excludesites", allparams)
append = parse_param("append", allparams)

num_refs = len(parse_ref_files("mixreference", allparams))
logfolder = parse_param("logfolder", allparams)



# def conRFMixAndMaskToBeagle(indfile_name, rephasedhaps_pref, em_iters, chroms):

concatenated = conRFMixAndMaskToBeagle(plinkbfile_pref + ".fam", plinkbfile_pref, em_iters, chroms)

# def filterForPCAmask(beaglefile, vitin, winin, fbkin, pcamask_refpanel, quality_cutoff, ancestry_number, ancestry_cutoff):

# Loop through all PCAmask reference panels that are specified in the param fil

pcamask_ref_panels = parse_ref_files("pcamaskrefpanel", allparams) # Ensure these are specified 
# return([outfile.name, vitout.name, winout.name, fbkout])
# in the same order as the normal reference panels
# def filterForPCAmask(beaglefile, vitin, winin, fbkin, pcamask_refpanel, quality_cutoff, ancestry_number, ancestry_cutoff):
pcamask_ready = []
for i,k in enumerate(pcamask_ref_panels):

	pcamask_ready.append(filterForPCAmask(concatenated[0], concatenated[1], concatenated[2], concatenated[3], num_refs, k, 
		qual_cut, i+1, anc_cut, window_size, exclude_inds, exclude_sites, append))
	
		
		
# Write temp file to signal master script
temp = open(plinkbfile_pref + "_cat", "w")
temp.write("All done.")
temp.close()
