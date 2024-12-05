#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-28

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/07_bighorn_rad

# define the reference genome
refgenome=/home/jmanthey/references/GCA_042477335.2_ARS-UI_OviCan_v2_genomic.fna 

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/*vcf.gz > \
${workdir}/04_vcf/${region_array}.vcf

# filter for 5 datasets

# dataset 1 for RAxML
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep dataset1.txt \
--max-missing 0.95 --mac 3 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/07_phylogeny/phylo_${region_array}

# dataset 2 for for observed heterozygosity and genetic distance
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep dataset2.txt \
--max-missing 0.95 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/08_stats/stats_${region_array}

# dataset 3 for pca of desert bighorn w/o Texas
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep dataset3.txt \
--max-missing 0.95 --mac 3 --max-alleles 2 --max-maf 0.49 --thin 10000 --recode \
--recode-INFO-all --out ${workdir}/06_structure/noTexas_${region_array}

# dataset 4 for pca and admixture of all desert bighorn
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep dataset4.txt \
--max-missing 0.95 --mac 3 --max-alleles 2 --max-maf 0.49 --thin 10000 --recode \
--recode-INFO-all --out ${workdir}/06_structure/desert_${region_array}

# dataset 5 for kinship and EEMS analyses (just Texas) 
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep dataset5.txt \
--max-missing 0.95 --mac 2 --max-alleles 2 --max-maf 0.49 --thin 10000 --recode \
--recode-INFO-all --out ${workdir}/09_relatedness/Texas_${region_array}

