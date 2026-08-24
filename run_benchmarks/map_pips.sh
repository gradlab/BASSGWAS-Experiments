mkdir $1;

for i in $(seq 0 20);
do
sh process_pips_R.sh $1 $2 $3 $i $i
sh process_pips_A.sh $1 $2 $3 $i $i
done
