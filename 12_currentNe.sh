# three steps


# Step 1. Interactive:

grep "Baylor" dataset5.txt > keep.baylor.txt
grep "Beaches" dataset5.txt > keep.beaches.txt
grep "BlackGap" dataset5.txt > keep.blackgap.txt
grep "ElephantMtn" dataset5.txt > keep.elephant.txt
grep "SierraDiablo" dataset5.txt > keep.diablo.txt
grep "VanHorn" dataset5.txt > keep.vanhorn.txt



# Step 2. Script submit:

#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=ne_filter
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-28

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/07_bighorn_rad

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# filter for 6 populations 
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keep.baylor.txt \
--max-missing 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --thin 100000 --recode \
--recode-INFO-all --out ${workdir}/10_popsize/baylor_${region_array}

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keep.beaches.txt \
--max-missing 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --thin 100000 --recode \
--recode-INFO-all --out ${workdir}/10_popsize/beaches_${region_array}

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keep.blackgap.txt \
--max-missing 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --thin 100000 --recode \
--recode-INFO-all --out ${workdir}/10_popsize/blackgap_${region_array}

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keep.elephant.txt \
--max-missing 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --thin 100000 --recode \
--recode-INFO-all --out ${workdir}/10_popsize/elephant_${region_array}

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keep.diablo.txt \
--max-missing 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --thin 100000 --recode \
--recode-INFO-all --out ${workdir}/10_popsize/diablo_${region_array}

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keep.vanhorn.txt \
--max-missing 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --thin 100000 --recode \
--recode-INFO-all --out ${workdir}/10_popsize/vanhorn_${region_array}



# Step 3. interactive:

# interactive session
interactive -p nocona -c 8

# set working directory
workdir=/lustre/scratch/jmanthey/07_bighorn_rad/10_popsize
cd $workdir

########################################
######## filtering
########################################

source activate bcftools

# remove sex chromosomes
rm *CM090947.1*
rm *CM090672.1*

# cat together
grep "^#" baylor_CM090671.1.recode.vcf > baylor.vcf

for i in $( ls baylor_*.recode.vcf ); do grep -v "^#" $i >> baylor.vcf; done

grep "^#" blackgap_CM090671.1.recode.vcf > blackgap.vcf

for i in $( ls blackgap_*.recode.vcf ); do grep -v "^#" $i >> blackgap.vcf; done

grep "^#" beaches_CM090671.1.recode.vcf > beaches.vcf

for i in $( ls beaches_*.recode.vcf ); do grep -v "^#" $i >> beaches.vcf; done

grep "^#" elephant_CM090671.1.recode.vcf > elephant.vcf

for i in $( ls elephant_*.recode.vcf ); do grep -v "^#" $i >> elephant.vcf; done

grep "^#" diablo_CM090671.1.recode.vcf > diablo.vcf

for i in $( ls diablo_*.recode.vcf ); do grep -v "^#" $i >> diablo.vcf; done

grep "^#" vanhorn_CM090671.1.recode.vcf > vanhorn.vcf

for i in $( ls vanhorn_*.recode.vcf ); do grep -v "^#" $i >> vanhorn.vcf; done


# remove the individual scaffold vcfs 
rm *recode*




# number of chromosomes in files (all had 26)
grep -v "#" baylor.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
wc -l chrom_map.txt  # 26
grep -v "#" beaches.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
wc -l chrom_map.txt  # 26
grep -v "#" blackgap.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
wc -l chrom_map.txt  # 26
grep -v "#" elephant.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
wc -l chrom_map.txt  # 26
grep -v "#" diablo.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
wc -l chrom_map.txt  # 26
grep -v "#" vanhorn.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
wc -l chrom_map.txt  # 26

# run the program
/home/jmanthey/currentNe/currentNe baylor.vcf 26

/home/jmanthey/currentNe/currentNe beaches.vcf 26

/home/jmanthey/currentNe/currentNe blackgap.vcf 26

/home/jmanthey/currentNe/currentNe elephant.vcf 26

/home/jmanthey/currentNe/currentNe diablo.vcf 26

/home/jmanthey/currentNe/currentNe vanhorn.vcf 26









