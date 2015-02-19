### Merging Beagle Files
### By Julian Homburger
import argparse
import sys


## Get both file names from command line
## Use options parse for -a, -b, -o 

parser = argparse.ArgumentParser()
parser.add_argument("-a", help="First Input File")
parser.add_argument("-b", help="Second Input File")
parser.add_argument("-o", help="Output file name")
# Arguments for Viterbi files, all three need to be specified for Viterbi combinations
parser.add_argument("-va", default="", help="Vit file for beagle A" )
parser.add_argument("-vb", default="", help="Vit file for beagle B" )
parser.add_argument("-vo", default="", help="Output Vit file")
parser.add_argument("-w", default='True', help="Should I write a windows file?")

args = parser.parse_args()

usevit = False
if (args.va != "" and args.vb != "" and args.vo != ""):
	usevit = True
	print "Combining Beagle and Viterbi files..."
else:
	print "Combining only Beagle files..."

bgl_one = open(args.a, "r")
bgl_two = open(args.b, "r")
bgl_out = open(args.o, "w")
## Read one file completely into hash, hash is:
bone_id = bgl_one.readline().strip()

if (usevit):
	vit_one = open(args.va, "r")
	vit_two = open(args.vb, "r")
	vit_out = open(args.vo, "w")

# Track the order of the rsIDs
bone_hash = dict()
i = 1
pos_hash_one = dict()
for bline in bgl_one:
	bsplits = bline.strip().split("\t", 2)
	bone_hash[bsplits[1]] = bsplits[2]
	pos_hash_one[bsplits[1]] = i
	i += 1
print i
# RS ID -> string of results

# Kind of memory intensive but meh

## Then, read the next file line by line
## If rsID is in the hash, then output this line to the output file
## If not, skip this line

btwo_id = bgl_two.readline().strip().split("\t",2)[2]

bgl_out.write(bone_id + "\t" + btwo_id + "\n")

# Now scan B2 to find matches
j = 1
pos_hash_two = dict()
final_order = list()
for b2l in bgl_two:
	b2s = b2l.strip().split("\t", 2)
	if b2s[1] in bone_hash:
		bgl_out.write(b2s[0] + "\t" + b2s[1] + "\t" + bone_hash[b2s[1]] + "\t" + b2s[2] + "\n")
		final_order.append(b2s[1])
	pos_hash_two[b2s[1]] = j
	j += 1
print j

## Now I have the position hash and the order of the rsIDs in the output file
## I also have the order of the individuals
if (not usevit):
	print "Finished combining beagle files, output to: " + args.o
	sys.exit()
else:
	print "Combining Viterbi files..."

## For each individual in viterbi one, read vit file and only keep positions relevant to the new beagle file

for v1 in vit_one:
	vsplit = v1.strip().split()
	vout =[vsplit[0]]
	for f in final_order:
		this_pos = pos_hash_one[f]
		#print this_pos
		#print len(vsplit)
		try:
			vout.append(vsplit[this_pos])
		except (IndexError):
			print this_pos
			print len(vsplit)
	vit_out.write("\t".join(map(str,vout)) + "\n")
	
	
	
## Do the same for the second viterbi file


for v2 in vit_two:
	vsplit = v2.strip().split()
	vout =[vsplit[0]]
	for f in final_order:
		this_pos = pos_hash_two[f]
		vout.append(vsplit[this_pos])
	vit_out.write("\t".join(map(str,vout)) + "\n")
	
### Generate a window file, default is true:

if (args.w == 'True'):
	wfile = open(args.vo + ".windows", "w")
	for k,f in enumerate(final_order):
		wfile.write("window" + str(k) + "\t" + f + "\n") 

	

