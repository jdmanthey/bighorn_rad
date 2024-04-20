#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/10_bighorn_rad

# run bcftools to merge the vcf files
bcftools merge -m id ${workdir}/03_vcf/*vcf.gz > ${workdir}/04_vcf/total.vcf

