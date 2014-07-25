# send this with qsub -v par="__",chr="___" 

if [ $# -eq 2 ]
then
	par=$1
	echo "singleChrRFMix.sh is taking arguments from command line input"
	chr=$2
fi
python singleChrRFMix.py $par $chr

