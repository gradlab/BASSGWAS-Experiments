#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=120G
#SBATCH -t 0-02:00
#SBATCH --job-name=annotate_pips
#SBATCH -o process_res.out
#SBATCH -e process_res.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 
#SBATCH -W

Rscript annotate/annotate_all.R $2 annotate/references.txt $1/annotated_hits.csv $1 
EOT
