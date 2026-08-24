DDIR=$1 #../data/lab_cip
PREFIX=$2
BATCH=$3
BEGIN=$4
END=$5
SEED=$6

PDIR=./BASSGWAS

export JULIA_NUM_THREADS=4
JULIA_NUM_THREADS=4

for i in $(seq $BEGIN $END);
do
    echo iteration $i
    IDIR=$PREFIX/RUN_$((i - 1))
    ODIR=$PREFIX/RUN_$i
    mkdir $ODIR

    julia $PDIR/probit_runner.jl --pfile $DDIR/patterns_centred.csv --ufile $DDIR/U.csv --sfile $DDIR/S.csv --obsfile $IDIR/data_bin.csv --odir $IDIR --batchsz $BATCH --seed $SEED$i --nthin 500 --nchains 4 --ndraws 1000
    cat $DDIR/data_bin.csv | (read line; echo "$line"; grep -wf $IDIR/design.txt) > $ODIR/data_bin.csv 
done
