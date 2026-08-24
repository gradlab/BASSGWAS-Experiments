#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -p shared,sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=60G
#SBATCH -t 0-03:00
#SBATCH --job-name=annotate_pips
#SBATCH -o annotate.out
#SBATCH -e annotate.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

rm remaining_kmers.fa
rm remaining_kmers.txt

Rscript annotate/annotate_all.R data/azi_mic_2/sequences.csv annotate/references.txt data/azi_mic_2/annotated_hits.csv data/azi_mic_2

rm remaining_kmers.fa
rm remaining_kmers.txt

Rscript annotate/annotate_all.R data/cip/sequences.csv annotate/references.txt data/cip/annotated_hits.csv data/cip

rm remaining_kmers.fa
rm remaining_kmers.txt

Rscript annotate/annotate_all.R data/tet10/sequences.csv annotate/references.txt data/tet10/annotated_hits.csv data/tet10

rm remaining_kmers.fa
rm remaining_kmers.txt

Rscript annotate/annotate_all.R data/tet2.5/sequences.csv annotate/references.txt data/tet2.5/annotated_hits.csv data/tet2.5

rm remaining_kmers.fa
rm remaining_kmers.txt
EOT
