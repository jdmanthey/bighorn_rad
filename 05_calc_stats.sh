#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=stats
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/07_bighorn_rad/08_stats

# define variables
vcf_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/vcf_stats.txt | tail -n1 )

# Rscript command for stats
Rscript _window_stat_calculations.r ${workdir}/${vcf_array} popmap_stats.txt


