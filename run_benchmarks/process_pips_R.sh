#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p hsph,shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=80G
#SBATCH -t 0-14:00
#SBATCH --job-name=process_res
#SBATCH -o process_res.out
#SBATCH -e process_res.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

parallel -j 6 "Rscript annotate/process_pips_2.R $2/{2}/RUN_{1} $3/annotated_hits.csv $3/sequences.csv $1/pips_{2}_{1}" ::: {$4..$5} ::: `find $2 -mindepth 1 -maxdepth 1 -type d -name R* -exec basename {} \; | xargs`;
EOT
