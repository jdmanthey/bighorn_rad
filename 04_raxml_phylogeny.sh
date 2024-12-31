#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=raxml
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=24
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

raxmlHPC-PTHREADS-SSE3 -T 24 -f a -x 50 -m ASC_GTRGAMMA --asc-corr=lewis -p 253 -N 100 \
-s /lustre/scratch/jmanthey/07_bighorn_rad/07_phylogeny/bighorn.fasta \
-n bighorn.tre \
-w /lustre/scratch/jmanthey/07_bighorn_rad/07_phylogeny/



