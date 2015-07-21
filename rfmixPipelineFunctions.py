### functions for calling RFMixPipeline

### This script holds all the functions used in the other RFMix pipeline scripts

### By J Homburger

import sys
import os
import re
import subprocess

## Function that takes a file and runs shapeIt - should return prefix of shapeit output files
## 

## TO DO:

## 1. Fix for multiple ancestries

def shapeItPhase(plink_name, genetic_map, known_haplotypes, max_threads=2, duohmm=False):
	out_pref = plink_name + "_shapeout"
	if (duohmm):
		thiscall = " ".join(["--input-bed " + plink_name + ".bed " + plink_name + ".bim " + plink_name + ".fam",
			"--input-map " + genetic_map, "--output-max " + out_pref, "-T " + str(max_threads), "--output-log " + plink_name + "_shapeitlog", "--duohmm"])
		subprocess.call("./shapeit " + thiscall, shell=True)
		#shape_log = subprocess.check_output(["shapeit" + thiscall])
		subprocess.call("rm " + plink_name + "_shapeitlog*", shell=True)	
	elif (known_haplotypes==""):
		# run command normally
		thiscall = " ".join(["--input-bed " + plink_name + ".bed " + plink_name + ".bim " + plink_name + ".fam",
			"--input-map " + genetic_map, "--output-max " + out_pref, "-T " + str(max_threads), "--output-log " + plink_name + "_shapeitlog"])
		print thiscall
		#os.system("./shapeit " + thiscall)
		subprocess.call("./shapeit " + thiscall, shell=True)
		#shape_log = subprocess.check_output(["shapeit" + thiscall])
		subprocess.call("rm " + plink_name + "_shapeitlog*", shell=True)	
	else:
		# run command with reference haplotypes
		ref_haps = known_haplotypes + ".haps"
		ref_samp = known_haplotypes + ".samp"
		ref_leg = known_haplotypes + ".leg"
		ref_call = "--input-ref " + " ".join([ref_haps, ref_samp, ref_leg])
		thiscall = " ".join(["--input-bed " + plink_name + ".bed " + plink_name + ".bim " + plink_name + " .fam",
			"--input-map " + genetic_map, "--output-max " + out_pref, "-T " + str(max_threads), ref_call])
		shape_log = subprocess.check_output(["shapeit", thiscall])
	return(out_pref)


#myout = shapeItPhase(input_plink, genetic_map, known_haplos_name)
#print myout
#ref_phased = []
#for anc in anc_refs:
#	ref_phased.append(shapeItPhase(anc, genetic_map, "")) # not using known haplotypes here

def readShapeItFile(shapeitfile):
	"""
	Helper function for below. This reads a shapeIt output
	file and returns two lists, one of the allele codes and one of the rsIDs/alleles
	"""
	alleles = []
	haps = []
	for sline in shapeitfile:
		splits = sline.strip().split()
		alleles.append([splits[0], splits[1], splits[2], splits[3], splits[4]])
		haps.append("".join(splits[5:]))
	return([haps, alleles])

def checkPolarity(a1, a2):
	"""
	Given allele vectors from the readShapeItFile function, checks if they are the same polarity
	"""
	assert(a1[1] == a2[1]) # Check to make sure these are the same
	if (a1[3] == a2[3] and a1[4] == a2[4]):
		return True
	elif (a1[3] == a2[4] and a1[4] == a2[3]):
		return False
	print a1
	print a2
	assert(False) # This means there is a strand issue	

def flip10(haps):
	"""
	Flips the 1/0 haplotype coding
	"""
	outhaps = ""
	for a in haps:
		if (a == "1"):
			outhaps += "0"
		elif (a == "0"):
			outhaps += "1"
	return(outhaps)

# this function takes as input both the admixed shapeit outputs and the ancestral shapeit outputs
def shapeItToRFMixMultiClass(rfmix_admix, rfmix_ancestrals):
	"""
	This function takes input in the form of a string for the admixed samples and a list of size > 1
	of the ancestral phased samples. Each input should be complete prefixes with filepaths for each
	"""
	print rfmix_admix + ".haps"
	hapfile = open(rfmix_admix + ".haps", "r")
	allele_file = open(rfmix_admix + ".alleles", "w") # File that stores 10 haplotypes
	mapfile = open(rfmix_admix + ".amaps", "w")
	# Need to keep track of population numbers here
	
	hap_strings, alleles = readShapeItFile(hapfile)
	pop_nums = [0]*len(hap_strings[0]) # Add admixed ind codes to pop_numbers array
	i = 1
	for thisanc in rfmix_ancestrals:
		thisfile = open(thisanc + ".haps", "r")
		this_haps, this_alleles = readShapeItFile(thisfile)
		pop_nums.extend([i]*(len(this_haps[0])))
		for j in range(0,len(hap_strings)): 
			# need to check if the alleles are flipped the same way
			if not checkPolarity(alleles[j], this_alleles[j]):
				this_haps[j] = flip10(this_haps[j])
			hap_strings[j] += this_haps[j]
		i += 1
	## Ok now I should have all the classes and haplotypes correctly
	classout = open(rfmix_admix + ".classes", "w")
	classout.write(" ".join(map(str,pop_nums)) + "\n")
	allele_file.write("\n".join(hap_strings) + "\n")
	for a in alleles:
		mapfile.write("\t".join(a) + "\n")
	allele_file.close()
	classout.close()
	mapfile.close()
	return([allele_file.name, classout.name, mapfile.name])

def makeSNPMap(snpfile, referencemap):
	"""
	Function that given a file of locus information and a reference genetic map.
	Creates a cM delimited mapfile for the given loci.
	Returns the name of the mapfile.
	"""
	bimfile = open(snpfile, "r") # open the input file
	mapfile = open(referencemap, "r")
	outfilename = re.sub(r'\.bim', '.markerpos', snpfile)
	posfilename = re.sub(r'\.bim', '.snp_locations', snpfile)
	outfile = open(outfilename, "w")
	posfile = open(posfilename, "w")
	# Initialize variables 
	previousCM = 0
	previousPos = 0
	i=0
	bimline = bimfile.readline().strip().split() # Pos 1 is rsID, Pos 3 is location
	for mapline in mapfile:
		if len(bimline) == 0:
			break		
		if i==0:
			i+=1
			continue
		mapline = mapline.strip().split()
		# Three cases: 1. SNP pos gt map pos
		while int(bimline[3]) < int(mapline[0]): # This means that the BIM file is behind the map file, so need to write output here with the interopolation
		# of the previous values
			diffCM = float(mapline[2]) - float(previousCM)
			diffpos = float(mapline[0]) - float(previousPos)
			multi = (float(bimline[3]) - float(previousPos))/diffpos
			cmout = multi*diffCM + float(previousCM)
			if cmout < 0: # this should not happen so if it does dump data and quit
				print i
				print cmout
				print diffCM
				print diffpos
				print previousCM
				print previousPos
				print bimline
				print mapline
				exit()

			outfile.write( str(cmout) +"\n")
			posfile.write( str(bimline[3]) + "\t" + str(cmout) + "\n")
			bimline = bimfile.readline().strip().split()
			if len(bimline) == 0:
				break		
		if len(bimline) ==0:
			break
		if bimline[3] == mapline[0]: # write out genetic position
			outfile.write( mapline[2]+ "\n")
			posfile.write( str(bimline[3]) + "\t" + mapline[2] + "\n")
			bimline = bimfile.readline().strip().split()
	
		#if bimline[3] > mapline[0]: # read next line in the map file
		#	previousCM = mapline[2]
		#	previousPos = mapline[0]
		#	continue
		# Hits this and continues if bimline is above mapline
		previousCM = mapline[2]
		previousPos = mapline[0]
		i += 1
	outfile.close()
	return(outfile.name)

## Function that runs RFMix to rephase alleles. Returns name of rephased alleles output
def runRFMixRephase(hapfile, classfile, snp_locations, window_size, generations, em_iters = 2, max_threads=2, trioPhase=False):
	"""
	This function takes as input the names of the concatenated haplotype file and the classfile
	It returns the name of the rephased allele file. 
	Note that here I do not include reference panels in the EM
	"""
	out_pref = re.sub(r".alleles", "", hapfile)
	os.chdir("RFMix") # RFMix doesn't like running outside its directory
	if (trioPhase):
		rfmix_call = " ".join(["python RunRFMix.py TrioPhased", hapfile,
			classfile, snp_locations, "--num-threads " + str(max_threads), 
			"-e " + str(em_iters), "-o " + out_pref, "--forward-backward"  ])
	else:
		rfmix_call = " ".join(["python RunRFMix.py PopPhased", hapfile,
			classfile, snp_locations, "--num-threads " + str(max_threads), 
			"-e " + str(em_iters), "-o " + out_pref, "--forward-backward"  ])
	print rfmix_call
	subprocess.call(rfmix_call, shell=True)
	os.chdir("../")
	return(out_pref)
	

def recodeAllele(allele, zero, ones):
	"""
	Recodes all 0's to the string zero and all 1's to the string 1
	If the allele is not a 0 or 1 - unlikely - it returns the input
	"""
	if allele=="0":
			return zero
	if allele=="1":
			return ones
	return allele

## Function that converts RFMix output to beagle output and transposes viterbi files.
def conRFMixAndMaskToBeagle(indfile_name, rephasedhaps_pref, em_iters, win_size, chroms):
	"""
	This takes the output of RFMix of the chromosomes in the chromosome vector and concatenates it
	into a beagle file, transposed viterbi file, and creates a new window file.
	Note that depending on the size of the RFMix file, this function may use a large amount of memory.
	"""
	### First get individual information
	window_id = 0
	em_iter = em_iters
	indfile = open(indfile_name, "r")	
	inds = []
	for line in indfile:
		splits = line.strip("\r\n").split()
		inds.append(splits[1] + "_A")
		inds.append(splits[1] + "_B")

	allloci = []
	outfilename = rephasedhaps_pref + "_w" + str(win_size) + ".beagle"
	outfile = open(outfilename, "w")
	outfile.write("I\tid\t" + "\t".join(inds) + "\n")
	## Write genotype data out to file

	vitout = open(rephasedhaps_pref + ".vit", "w")
	winout = open(rephasedhaps_pref + ".windows", "w")
	fbkout = rephasedhaps_pref + ".fbk"
	if os.path.exists(fbkout):
		os.remove(fbkout)
	vitlist = []
	for chrom in chroms:
		print chrom
		shapeitfilename = rephasedhaps_pref + "_chr" + str(chrom) + "_shapeout.allelesRephased" + str(em_iters) + ".txt"
		shapeitfile = open(shapeitfilename, "rb")
		fbkin_name = rephasedhaps_pref + "_chr" + str(chrom) + "_shapeout." + str(em_iters) + ".ForwardBackward.txt"
		os.system('cat ' + fbkin_name + " >> " + fbkout) # Concatenate files together
		markerin = rephasedhaps_pref + "_chr" + str(chrom) + "_shapeout.amaps"
		markerfile = open(markerin, "r")
		loci=[]
		alleles = {}
		for mline in markerfile:
			msplit = mline.strip().split()
			loci.append(msplit[1])
			alleles[msplit[1]] = [msplit[3], msplit[4] ]

		allloci.extend(loci)
		for j,line in enumerate(shapeitfile):
			sline = line.strip("\r\n")
			zero, ones = alleles[loci[j]]
			fixed = [ recodeAllele(k, zero, ones) for k in sline ]
			outfile.write("M\t" + loci[j] + "\t" + "\t".join(fixed) + "\n")
		vitfile = open(rephasedhaps_pref + "_chr" + str(chrom) + "_shapeout." + str(em_iters) + ".Viterbi.txt", "r")
		vitlist.extend([x.strip().split() for x in vitfile])
		shapeitfile.close()
		vitfile.close()
		
	# This will transpose the whole Viterbi file
	# Yikes this may take a lot of memory
	for i,x in enumerate(zip(*vitlist)):
		vitout.write(inds[i] + "\t")
		for y in x:
			vitout.write(y+"\t")
		vitout.write("\n")
		### This doesn't quite work yet so make sure to fix it next time
	for l in allloci:
		winout.write("window" + str(window_id) + "\t" + l + "\n")
		window_id += 1
	return([outfile.name, vitout.name, winout.name, fbkout])


def singcomp(char):
	'''
	Given a single DNA base, returns its complement
	'''
	if char == "A":
		return("T")
	elif char == "T":
		return("A")
	elif char == "C":
		return("G")
	elif char == "G":
		return("C")
	elif char == "1" or char == "0": # This function won't change 1/0 codings
		return char
	print char
	raise ValueError("Not a recognized DNA nucleotide")
	
	
# Essentially copy your PCAdmix filtering script for this - can reuse a lot of that code
def filterForPCAmask(beaglefile, vitin, winin, fbkin, num_refs, pcamask_refpanel, quality_cutoff, ancestry_number, ancestry_cutoff, window_size, exclude_inds = "", exclude_sites = "", append=""):
	'''
	This function takes input in the form of concatenated beagle, viterbi, window, and forward-backwards files
	and creates inputs to PCAMask that have markers that overlap with the ancestral reference panels.
	This can be a different reference panel than was used for local ancestry calling earlier. 
	'''
	missing_mask = "-9"
	beaglefile = beaglefile # Name of the admixed beagle file
	vitfilename= vitin # Name of the admixed viterbi calls
	ancbeagle = pcamask_refpanel # Name of the ancestral beagle file
	fbkname = fbkin # Name of the forward backward probabilities for the admixed viterbi calls
	windowname= winin # Name of the window file produced by PCAdmix
	prefix = re.sub(r".beagle", "", beaglefile)
	vitoutname= prefix + "q" + str(int(quality_cutoff*100)) + "_b" + str(int(100*ancestry_cutoff)) + "_w" + str(int(100*window_size)) + "_a" + str(ancestry_number) + append + ".vit" #Name of output viterbi file
	winoutname=prefix + "q" + str(int(quality_cutoff*100)) + "_b" + str(int(100*ancestry_cutoff)) + "_w" + str(int(100*window_size)) + "_a" + str(ancestry_number) + append + ".windows" # Name of output window file
	outname = prefix + "q" + str(int(quality_cutoff*100)) + "_b" + str(int(100*ancestry_cutoff)) + "_w" + str(int(100*window_size)) + "_a" + str(ancestry_number) + append + ".beagle" # Name of output beagle file
	ancmarkers= {} # Holds all the ancestral markers and their alleles
	j = 0 
	target_pop = str(ancestry_number)
	exclude = set()
	if exclude_sites != "":
		exlude_sites = open(exclude_sites, "r")
		for e in exclude_sites:
			exclude.add(e.strip())
	
	inds_to_exclude = set()
	if exclude_inds != "":
		exclude_inds = open(exclude_inds, "r")
		for y in exclude_inds:
			inds_to_exclude.add(y.strip())
	
	block_cutoff = ancestry_cutoff
	ancbeagle = open(ancbeagle, "r")
	beaglefile = open(beaglefile, "r")
	fbkcut = float(quality_cutoff)
	
	for ancline in ancbeagle:
		if j == 0:
			j+=1
			continue
		ancsplit = ancline.strip().split()
		if ancsplit[1] in exclude:
			continue

		myset = set(ancsplit[2:])
		if len(myset) != 2: # remove markers monomorphic in reference panel
			continue
		ancmarkers[ancsplit[1]] = myset

	outfile = open( outname, 'w')
	in_windows=set()
	windowfile= open(windowname, "r")
	winout = open(winoutname, "w")
	windows = []
	block_size=[]

	for winline in windowfile:
		winsplit = winline.strip().split()
		keep = [winsplit[0]]
		for w in winsplit[1:]:
			if w in ancmarkers:
				keep.append(w)
				in_windows.add(w)
		if len(keep) > 1:
			winout.write("\t".join(keep) + "\n")
			windows.append(1)
			block_size.append(len(keep)-1)
		else:
			windows.append(0)
			block_size.append(0)

	in_anc = set(ancmarkers.keys())
	both = in_anc & in_windows
	index_to_snp = []
	i=0
	k=0
	outfile_lines = []
	print ("Reading beagle file into memory for block check")
	for line in beaglefile:
		if i == 0:
			i+=1
			line = line.strip().split()
			outfile_lines.append(line)
			continue
		line = line.strip().split()
	
		# map snp index back to rsid - will check if it's in ancmarkers later
		index_to_snp.append(line[1])
		if line[1] in both:
			thisset = set(line[2:])
		
			if not ancmarkers[line[1]] == thisset:
				if len(thisset) == 2:
					line[2:] = [singcomp(a) for a in line[2:]]
					k += 1
				if len(thisset) != 2:
					thiscomp = singcomp(thisset.pop())
					if thiscomp in ancmarkers[line[1]]:
						line[2:] = [singcomp(a) for a in line[2:]]
				if not set(line[2:]) <=	ancmarkers[line[1]]:
					print set(line[2:])
					print ancmarkers[line[1]]
					sys.exit(1)
			outfile_lines.append(line)
		i+=1

	print ( "I found and corrected " + str(k) + " strand issues.")

	print("Writing filtered viterbi file...")
	vitout = open(vitoutname, "w")
	viterbifile = open(vitfilename, "r")
	fbk = open(fbkname, "r")

	inds_to_keep = [1,1] # keep first two columns

	miss = 0
	totalblocks=0
	
	probs = []
	### fix this part 
	### get matrix of fbk probabilities
	### matrix has columns of individuals and rows of alleles
	for y, fline in enumerate(fbk): # for each row
		fline1 = fline.strip().split()
		probs.append([]) # append an extra array for each row
		for r in range(0,len(fline1), num_refs):
			probs[y].append(max(map(float,fline1[r:(r+num_refs)])))
	
	### Viterbi file has rows of individuals and columns of alleles
	for g, vline in enumerate(viterbifile):
		vline = vline.strip().split()
		fline1 = fbk.readline().strip().split() # read fbk file in line with Viterbi file

		id = vline[0]
		anccalls = vline[1:]
		output_calls = []

		for i in range(0,len(windows)): 
			if windows[i] == 1:
				if float(probs[i][g]) < fbkcut: # check to see if it passes quality theshold
					output_calls.append(missing_mask)
					miss +=1 
					totalblocks+=1
				else:
					output_calls.append(anccalls[i])
					totalblocks+=1
		
		counts = output_calls.count(target_pop)
		if (len(output_calls) < 1):
			inds_to_keep.append(0)
			continue
		if (float(counts)/float(len(output_calls))) < block_cutoff or id in inds_to_exclude:
			inds_to_keep.append(0)
		else: # write to output file and save position of individual to keep
			inds_to_keep.append(1)
			vitout.write(id + "\t" + "\t".join(map(str,output_calls)) + "\n")

			
		#vitout.write('%s\t%s\n' % (id, '\t'.join([anccalls[i] for i in range(len(anccalls)) if index_to_snp[i] in ancmarkers])))

	print("Performing block check and outputting beagle file")
	sum = 0
	for d in inds_to_keep:
		sum += d
	inds_removed = len(inds_to_keep) - sum
	print ( str(inds_removed) + " individuals of " + str(len(inds_to_keep) -2) + " were removed for low numbers of ancestral regions")
	print ( str(miss) + " blocks out of " + str(totalblocks) + " were masked due to low quality.")

	flag=0
	for thisout in outfile_lines:
	
		keep = [ thisout[m] for m in range(len(thisout)) if inds_to_keep[m] == 1 ]
		if flag ==0:
			flag=1
		outfile.write("\t".join(keep) + "\n")

	return([outfile.name, vitout.name, winout.name])

def parse_param(paramstring, thisfile):
	'''
	Returns string of the parameter
	'''
	try:
		this = re.search(paramstring + r":.+\n", thisfile).group(0)
		this = re.sub(paramstring+r":", "", this).strip()
	except AttributeError:
		return ""
	return this

def parse_ref_files(paramstring, thisfile):
	'''
	Special parser for names of references because there can be multiple references
	'''
	try:
		refs = []
		pattern = paramstring + r":.+\n"
		these = re.findall(pattern, thisfile)
		for r in these:
			refs.append(re.sub(paramstring + r":", "", r).strip())
		return(refs)
	except AttributeError:
		return("")
	
