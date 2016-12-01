DIR=$1
RESULT_FILES=`ls $DIR*.results`

for file in $RESULT_FILES;
do

echo "python grapher.py $file $'\t'"

done
