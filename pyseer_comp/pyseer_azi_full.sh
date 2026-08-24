#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH -t 0-14:00
#SBATCH --job-name=pyseer_cip
#SBATCH -o pyseer_cip.out
#SBATCH -e pyseer_cip.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

bash ./process_full.sh azi azi_mic_2 
EOT
