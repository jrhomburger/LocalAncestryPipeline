### Attach Viterbi calls to IBD files

## Inputs:

# 1. Beagle file - to get locations of each RSID
# 2. Viterbi file for each individual with the local ancestry calls
# 3. Match file detailing beginning and end of the matches

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Germline Match Output File")
parser.add_argument("-m", help="PLINK map file with markers ordered")
parser.add_argument("-v", help="Viterbi output File from RFMix")
## Note outputting all the data will make a very large file
parser.add_argument("--output_string", help="Should I output the entire Viterbi string in the match file?" )

### First, get positions of each marker in a dictionary for easy lookup:
### Ie., read and store the order of the map file
args = parser.parse_args()

# Open all the files (crash here if you can't find one)
mapfile = open(args.m, "r")
germfile = open(args.g, "r")
vitfile = open(args.v, "r")
outfile = open(args.g + ".append", "w")

## Read markers file for positions

print "Reading map file..."

m_hash = dict()

for i,mline in enumerate(mapfile):
	msplit = mline.strip().split()
	m_hash[msplit[1]] = i
	
## Read in the entire Viterbi file for easy access
## Ah I love having lots of RAM

vhash = dict()

print "Reading Viterbi file..."

for vline in vitfile:
	vsplit = vline.strip().split("\t")
	vhash[vsplit[0]] = vsplit[1:]


def vitCounter(vitvec):
	"""
	Given a viterbi vector, counts the numbers of 1,2,3, and returns a list
	"""
	num1 = vitvec.count("1")
	num2 = vitvec.count("2")
	num3 = vitvec.count("3")
	return([num1, num2, num3])
	
def discordantVec(vec1, vec2):
	assert len(vec1) == len(vec2)
	discs =0
	for i in range(len(vec1)):
		if (vec1[i] != vec2[i]):
			discs+=1
	return(discs)

def mostCommonAnc(vec, anc):
	index = vec.index(max(vec))
	return(anc[index])

## Ok, now the fun/tricky part
## Read through the match file line by line
## Then, identify the start and stop for each line
## Using the position data, use the start and stop to pull the correct Viterbi regions
## Concatenate the viterbi calls into a comma separated string and summary stats
## Output string and summary stats to file

print "Processing match file..."

header = ["FamID1", "SampleID1", "FamID2", "SampleID2", "Chr", "Match_Start", "Match_End",
	"RS_start", "RS_end", "Number_SNPs", "GenDist", "Units", "Mismatches" , "homozygous1", "homozygous2",
	"Eur1", "Nat1", "Afr1", "Eur2", "Nat2", "Afr2", "Discordant", "Ancestry1", "Ancestry2"]
ident_start=0
outfile.write("\t".join(header) + "\n")
for line in germfile:
	splits = line.strip().split()
	## 0 FamID, SampID, FamID, SampID, Chr, # Pos start, Pos end, rs start, rs end, # snps, 
	
	sampID_1 = splits[1].replace(".1", "_B").replace(".0", "_A")
	sampID_2 = splits[3].replace(".1", "_B").replace(".0", "_A")

	rsStartPos = m_hash[splits[7]]
	rsEndPos = m_hash[splits[8]]
	
	if (rsStartPos == rsEndPos):
		ident_start += 1
		continue
	
	### Get the two vectors with the local ancestry calls:
	
	vitv1 = vhash[sampID_1][rsStartPos:(rsEndPos+1)]
	vitv2 = vhash[sampID_2][rsStartPos:(rsEndPos+1)]
	ancs = ["European", "American", "African"]
	
	### Compare them:
	### I want, 1. Count Eur Calls, 2. Count Afr Calls, 3. Count Nat Calls for each individual
	
	counts1 = vitCounter(vitv1)
	counts2 = vitCounter(vitv2)
	
	max1 = mostCommonAnc(counts1, ancs)
	max2 = mostCommonAnc(counts2, ancs)

	
	## Also want the number of discordant local ancestry calls
	discs = discordantVec(vitv1, vitv2)
	toAppend = "\t".join(map(str, counts1)) + "\t" + "\t".join(map(str, counts2)) + "\t" + str(discs) + '\t' + max1 + '\t' + max2

	if (args.output_string):
		toAppend = toAppend + "\t" + ",".join(vitv1) + "\t" + ",".join(vitv2)
	
	writeoutput = "\t".join(splits) + "\t" + toAppend + "\n"
	outfile.write(writeoutput)
	

## Ok I think that should do it, then we can compare the European calls and the lengths of the IBD matches
print ("I found " + str(ident_start) + " positions with identical starts")



	
	
	


