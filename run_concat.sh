module load python/2.7

if [ -e /srv/gs1/software/python/2.7.6/bin/python ]
then
	alias python='/srv/gs1/software/python/2.7.6/bin/python'
fi
if [ $# -eq 1 ]
then
	par=$1
fi


if [ -e /srv/gsfs0/software/python/2.7.6/bin/python ]                
then                                                                                                                                                                                                      
        /srv/gs1/software/python/2.7.6/bin/python concatRFMixOut.py $par 
else
	python concatRFMixOut.py $par                                                                                                                                           
fi  
 
