#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=60G
#SBATCH -t 3-00:00
#SBATCH --job-name=benchmark.adapt.$2
#SBATCH -o bench.adapt.$2.out
#SBATCH -e bench.adapt.$2.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

sh scripts/run_adaptive_old.sh $1 $2 $3 $4 $5 $6
EOT
