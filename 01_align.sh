#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-244

module load intel java bwa singularity

export SINGULARITY_CACHEDIR="/lustre/work/jmanthey/singularity-cachedir"

# define main working directory
workdir=/lustre/scratch/jmanthey/10_bighorn_rad

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome
refgenome=/home/jmanthey/references/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna

# run bbduk to clean and filter
/lustre/work/jmanthey/bbmap/bbduk.sh in1=${workdir}/00_fastq/${basename_array}_R1.fastq.gz in2=${workdir}/00_fastq/${basename_array}_R2.fastq.gz out1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz out2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz minlen=50 ftl=5 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bwa mem
bwa mem -t 4 ${refgenome} ${workdir}/01_cleaned/${basename_array}_R1.fastq.gz ${workdir}/01_cleaned/${basename_array}_R2.fastq.gz > ${workdir}/01_bam_files/${basename_array}.sam

# convert sam to bam
samtools view -b -S -o ${workdir}/01_bam_files/${basename_array}.bam ${workdir}/01_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01_bam_files/${basename_array}.sam

# clean up the bam file
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk CleanSam -I ${workdir}/01_bam_files/${basename_array}.bam -O ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk SortSam -I ${workdir}/01_bam_files/${basename_array}_cleaned.bam -O ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam --SORT_ORDER coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk AddOrReplaceReadGroups -I ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam -O ${workdir}/01_bam_files/${basename_array}_final.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${basename_array}_final.bam
