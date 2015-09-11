module load python/2.7

if [ $# -eq 1 ]
then
	par=$1
fi 
python requality.py $par
