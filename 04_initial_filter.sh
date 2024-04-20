#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/10_bighorn_rad

# filter based on missing data (max = 10 individuals or 0.042 proportion)
vcftools --vcf ${workdir}/04_vcf/total.vcf --max-missing 0.958 \
--max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/total

