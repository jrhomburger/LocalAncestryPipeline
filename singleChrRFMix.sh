# send this with qsub -v par="__",chr="___" 
if [ -e /srv/gs1/software/python/2.7.6/bin/python]
then
        alias python='/srv/gs1/software/python/2.7.6/bin/python'
fi
if [ $# -eq 2 ]
then
	par=$1
	echo "singleChrRFMix.sh is taking arguments from command line input"
	chr=$2
fi
python singleChrRFMix.py $par $chr

