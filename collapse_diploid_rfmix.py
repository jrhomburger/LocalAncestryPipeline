### Collapse RFMix files into diploid bed file of ancestries

from optparse import  OptionParser

### Need as input
# 1. FBK file
# 2. 


def splitstr(option, opt, value, parser):
  return(setattr(parser.values, option.dest, value.split(',')))

USAGE = """
collapse_ancestry.py    --rfmix
                        --snp_locations
                        --ref_run
                        --sa_run
                        --vit
                        --fbk
                        --ind
                        --ind_info
                        --pop_labels
                        --out
"""

parser = OptionParser(USAGE)

parser.add_option('--rfmix', default='/home/armartin/rare/chip_collab/AFR_EUR_NATAM_') #ACB_chr1.rfmix.0.Viterbi.txt
parser.add_option('--snp_locations', default='/home/armartin/rare/chip_collab/scripts/chr')
parser.add_option('--vit', default='0')
parser.add_option('--fbk', default=None)
parser.add_option('--admixed_pop', default='ACB')
parser.add_option('--ind_info') #Individual IDs listed in the order they appear in RFMix output
parser.add_option('--pop_labels', type='string', action='callback', callback=splitstr, default=['AFR','EUR','NAT'],
                  help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
parser.add_option('--out', default='/home/armartin/rare/chip_collab/lai_plots/')
parser.add_option('--fbkcut', default='0.9', type=float)

(options,args)=parser.parse_args()

### First, find the individuals location

### Then, extract their genotypes

### Extract forward backward probabilities

### Identify switchpointd

### Write bedfile

chrs = range(1,23)
print options.pop_labels
pop_labels = options.pop_labels
fbk_threshold = options.fbkcut

#load parameters and files
admixed_pop = options.admixed_pop
vit = options.vit

ind_info = open(options.ind_info)
ind_list = []
for line in ind_info:
    myLine = line.strip().split()
    ind_list.append(myLine[0])

print ind_list

if options.ind is None:
    current_ind = ind_list
else:
    current_ind = options.ind

#note, incorporate admixture info here
#anc_info = open('/home/armartin/sa_analysis/data/covariates/pheno_ancestry.txt')
#anc_dict = {}
#anc_info.readline()
#for line in anc_info:
#    myLine = line.strip().split()
#    anc_dict[myLine[0]] = [myLine[-2], myLine[-6], myLine[-4], myLine[-5]]    

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def collapse_bed(current_ind):
	
	bedfile_out = open(options.out + current_ind + '.bed', 'w')
	    counter = 0
    for chr in chrs:
        print chr
        rfmix = open(options.rfmix + '_chr' + str(chr) + '_shapeout.' + vit + '.Viterbi.txt') #1.rfmix.0.Viterbi.txt
        if options.fbk is not None: #will need to change global proportion calculator, make new anc possibility (-9?), and color black
            fbk = open(options.rfmix +  '_chr' + str(chr) + '_shapeout.' + vit + '.ForwardBackward.txt')
        snp_locations = open(options.snp_locations + '_chr' + str(chr) + '.markerpos')
        snp_map = open(options.snp_locations + '_chr'+ str(chr) + '.snp_locations') #map of physical position -> genetic position
        ind = ind_list.index(current_ind)
        
               last_anc = None
        last_plot_bound = None
        last_hapa_anc = 0
        last_hapa_pos = None
        last_hapb_anc = 0
        last_hapb_pos = None
        last_switch = 0
        
        for line in rfmix:
            myLine = line.strip().split()
            my_pos = snp_locations.readline().strip()
            my_map = snp_map.readline().strip().split()
            #finds max of 3 ancestry posterior probabilities
            if options.fbk is not None:
                fbk_line = fbk.readline().strip().split()
                fbk_line = map(float, fbk_line)
                fbk_max = []
                for i in grouper(len(pop_labels), fbk_line):
                    fbk_max.append(max(i))
            #print fbk_max
            current_anc = [myLine[ind*2], myLine[ind*2+1]]
            current_info = [my_pos, myLine[ind*2], myLine[ind*2+1]]
            if fbk_max[ind*2] >= fbk_threshold:
                current_hapa_anc = myLine[ind*2]
            else:
                current_hapa_anc = -9
                current_anc[0] = -9
            if fbk_max[ind*2+1] >= fbk_threshold:
                current_hapb_anc = myLine[ind*2+1]
            else:
                current_hapb_anc = -9
                current_anc[1] = -9
            #current_hapb_anc = myLine[ind*2+1]
            current_hapa_pos = my_map[0]
            current_hapb_pos = my_map[0]
            current_hapa_cm = my_map[1]
            current_hapb_cm = my_map[1]
            current_plot_bound = my_pos
            if last_hapa_anc == 0:
                last_hapa_anc = current_hapa_anc
                last_hapa_pos = current_hapa_pos
                last_hapa_cm = current_hapa_cm
            if last_hapb_anc == 0:
                last_hapb_anc = current_hapb_anc
                last_hapb_pos = current_hapb_pos
                last_hapb_cm = current_hapb_cm
            if current_anc == last_anc:
                pass
            else:
                print [last_anc, chr, last_plot_bound, current_plot_bound, current_info, [current_hapa_anc, current_hapb_anc]]
                #print chr
                #print last_anc
                if last_plot_bound is not None:
                   # lai_proportions = track_lai_proportions(last_plot_bound, current_plot_bound, [current_plot_bound, last_anc[0], last_anc[1]], lai_proportions)
                    counter = counter + 2*(float(current_plot_bound) - float(last_plot_bound))
                    
                    
                    ### this is the meat here
                    ### If you find a switchpoint
                    if (last_hapa_anc is not None and last_hapb_anc is not None) and (last_hapa_anc != current_hapa_anc 
                    	or !last_hapb_anc != current_hapb_anc):
                    	###
                        if (last_hapa_anc != -9 and last_hapb_anc != -9) and last_hapa_anc != current_hapa_anc:
                            bedfile_out.write(str(chr) + '\t' + last_switch + '\t' + current_hapa_pos + '\t' + pop_labels[int(last_hapa_anc)-1] + ":" + pop_labels[int(last_hapb_anc)-1] + '\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
                       		last_switch = current_hapa_pos
                       		last_hapa_anc = current_hapa_anc
                        	last_hapa_pos = current_hapa_pos
                        else if (last_hapa_anc != -9 and last_hapb_anc != -9) and last_hapb_anc != current_hapb_anc:
                        	bedfile_out.write(str(chr) + '\t' + last_switch + '\t' + current_hapb_pos + '\t' + pop_labels[int(last_hapa_anc)-1] + ":" + pop_labels[int(last_hapb_anc)-1] + '\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
							last_switch = current_hapb_pos
							last_hapb_anc = current_hapb_anc
							last_hapb_pos = current_hapb_pos
                        else if (last_hapa_anc != -9) and last_hapa_anc != current_hapa_anc:
                        	bedfile_out.write(str(chr) + '\t' + last_switch + '\t' + current_hapa_pos + '\t' + pop_labels[int(last_hapa_anc)-1] + ":" + "UNK" + '\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
                       		last_switch = current_hapa_pos
						else if (last_hapb_anc != -9) and last_hapb_anc != current_hapb_anc:
                        	bedfile_out.write(str(chr) + '\t' + last_switch + '\t' + current_hapb_pos + '\t' + "UNK" + ":" + pop_labels[int(last_hapb_anc)-1] + '\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
							last_switch = current_hapb_pos
                        else:
                            bedfile_out.write(str(chr) + '\t' + last_hapa_pos + '\t' + current_hapa_pos + '\tUNK:UNK\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
                        
                last_plot_bound = current_plot_bound
                last_anc = current_anc #might need to change this for plotting purposes
        if current_hapa_anc == -9:
			a_anc = "UNK"
        else:
        	a_anc = pop_labels[int(current_hapa_anc)-1]
        if current_hapb_anc == -9:
        	b_anc = "UNK"
            hap_b.write(str(chr) + '\t' + last_hapb_pos + '\t' + current_hapb_pos + '\tUNK\t' + last_hapb_cm + '\t' + current_hapb_cm + '\n')
        else:
        	b_anc = pop_labels[int(current_hapb_anc)-1]
		bedfile_out.write(str(chr) + '\t' + last_hapb_pos + '\t' + current_hapb_pos + '\t' +  + '\t' + last_hapb_cm + '\t' + current_hapb_cm + '\n')
        print [last_anc, chr, last_plot_bound, current_plot_bound, current_info]
        #plot_chromosomes(last_anc, chr, last_plot_bound, current_plot_bound, current_info, ax)
    
    bedfile_out.close()
    
    

ind_info = open(options.ind_info)
ind_list = []
for line in ind_info:
    myLine = line.strip().split()
    ind_list.append(myLine[0])

print ind_list


for ind in ind_list:
	collapse_bed(ind)

    
	
	
