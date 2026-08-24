BEGIN=$1
END=$2
DIR=$3
DDIR=$4
NREP=$5
BATCH=$6

for i in $(seq 1 $NREP);
do
    sh scripts/run_rep_adaptive_old.sh $DDIR $DIR/A$i $BATCH $BEGIN $END 2$i$i$i
done
