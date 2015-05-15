if [ -e /srv/gs1/software/python/2.7.6/bin/python]
then
	alias python='/srv/gs1/software/python/2.7.6/bin/python'
fi
if [ $# -eq 1 ]
then
	par=$1
fi 
python concatRFMixOut.py $par
