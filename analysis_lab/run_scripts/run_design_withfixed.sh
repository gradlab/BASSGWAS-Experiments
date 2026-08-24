#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p hsph,shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --mem=120G
#SBATCH -t 1-00:00
#SBATCH --job-name=design1_long
#SBATCH -o design_l.out
#SBATCH -e design_l.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu

export JULIA_NUM_THREADS=4

julia ../BASSGWAS_oldprobit_runner.jl  --pfile data_with_bad/patterns_centred.csv --ufile data_with_bad/U.csv --sfile data_with_bad/S.csv --obsfile $1 --odir $2 --batchsz 4 --nthin 750 --nchains 8 --ndraws 1500 --inclnext $3 --seed $4
EOT
