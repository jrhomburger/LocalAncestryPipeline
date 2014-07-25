# README for Local Ancestry Pipeline #

By Julian Homburger
jhomburg@stanford.edu

### Local Ancestry Pipeline ###

This code is used to automate local ancestry calling. It starts with a set of PLINK bed files
and outputs local ancestry calls that are ready to enter PCAmask.

version 0.1


### How to use these scripts ###

There are multiple scripts that are controlled by a paramfile, much as in EIGENSOFT. 
To use, call python fullRFMixMaster.py paramfile.
The functions are divided into multiple scripts to facilitate parallelization
and their use on an SGE cluster.

REQUIREMENT:
The scripts expect the RFMix folder to be located in the same folder as the script.
The RFMix folder's name should be changed from RFMix_v1.x.x to RFMix.
If this is not done, the scripts will be unable to call RFMix.
The script also requires the shapeit executable file to be in the script directory.

#### So make sure you have in the same directory: ####
1.	Folder named RFMix with (compiled) RFMix programs in it
2.	ShapeIt executable

### Summary of components ###

* fullRFMixMaster.py - master script, called by the user and coordinates the other scripts
* rfmixPipelineFunctions.py - script that contains the functions used by the called scripts
* singleChrRFMix.sh - shell wrapper for singleChrRFMix.py
* singleChrRFMix.py - takes data from a single chromosome and runs it through ShapeIt and RFMix
* run_concat.sh - shell wrapper for concatRFMixOut.py
* concatRFMixOut.py - script that concatenates output into beagle files and creates files ready for ASPCA
* requality.sh - shell wrapper for requality.py 
* requality.py - script the runs when changing quality or ancestry cutoffs

### Parameter Descriptions and Codes ###

*** Note that it is imperative the parameter names in the paramfile
match those below exactly or the program will ignore the parameter ***

* plinkbfile: prefix of the plink bed file of admixed individuals (for instance, for 'lion.bed', use 'lion')
* mixreference: prefix of the plink bed file of continental reference panel. specify this once for each continental panel.
* pcamaskrefpanel: beagle file of subcontinental reference panel for PCAmask. Can be specified multiple times.
*** Note that currently the script assumes that the mixreference and pcamaskrefpanel parameters are in the same order
So, if your European ref panel is first on mixreference than you subcontinental European panel should be the first
specified on the pcamaskrefpanel. Suggestions for creative ways to link the two together are appreciated.***
* geneticmap: location of the geneticmap files. Replace the individual chromosome number with NNN. for example, if
we have the first chromosome map specified as path/geneticmapchr1b36.txt, then specify path/geneticmapchrNNNb36.txt
* maxthreads: max number of threads for each instance of singleChrRFMix.py. Default is 1.
* qsub: if 'True', it will use qsub. If anything else or omitted, will run locally.
* windowsize: RfMix window size, in centimorgans, default is 0.2
* generations: RfMix generations, default is 8
* emiterations: RfMix EM iterations, default is 2
* qualitycutoff: cutoff for quality when parsing for PCAmask - default is 0.95
* ancestrycutoff: cutoff for removing individuals with low ancestral pop ancestry, default is 0.25
* phasemix: If set to False, will not run phasing/RFMix portion of script. Default is True. 
* runconcat: If set to False, the concatenation and filtering steps will not be run
* requality: If set to True, neither the concatenation nor the phase/RFMix part of the script will run. 
Instead, it will recreate the PCAmask input files with the specified quality parameters
* chrom: Give a list of chr x,y,z (no spaces) will run only those chromosome numbers
* exclude_inds: file with each line a sampleID. These sample IDs will be excluded from PCAMask output
* exclude_sites: file with each line a locus ID. These loci will be excluded from PCAMask output
* append: appends an extension on the filtered PCAmask output files. Useful for telling different quality
filters apart. 
* skip_anc: use in conjunction with requality == TRUE. If inputted, will skip recreating the PCAmask files for
the given ancestry numbers. Note that you still will need the pcamaskrefpanel parameter specified for each ancestry -
however skipped ancestries can be given any file path
* logfolder: specifies a path for log files 



