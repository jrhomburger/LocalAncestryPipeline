module load python/2.7

# send this with qsub -v par="__",chr="___" 
if [ -e /srv/gsfs0/software/python/2.7.6/bin/python ]
then
        alias python276='/srv/gsfs0/software/python/2.7.6/bin/python'
fi
if [ $# -eq 2 ]
then
	par=$1
	echo "singleChrRFMix.sh is taking arguments from command line input"
	chr=$2
fi

if [ -e /srv/gsfs0/software/python/2.7.6/bin/python ]                                                                                                                                                      
then    
	echo "calling specific python instance" 
                                                                                                                                                          
	/srv/gsfs0/software/python/2.7.6/bin/python singleChrRFMix.py $par $chr  
else
	python singleChrRFMix.py $par $chr
fi 


