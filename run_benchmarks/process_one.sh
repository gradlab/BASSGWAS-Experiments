#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p hsph,shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH -t 0-4:00
#SBATCH --job-name=process_one
#SBATCH -o process_one.out
#SBATCH -e process_one.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

mkdir $1
Rscript annotate/process_pips_2.R $2 $4/annotated_hits.csv $3 $1/pips
EOT
