# interactive session with 4 cpu and 32 GB memory

interactive -p nocona -c 4 -m 8G

source activate bcftools

workdir=/lustre/scratch/jmanthey/07_bighorn_rad/

cd ${workdir}/06_structure

########################################
######## combine files
########################################

# remove sex chromosomes
rm noTexas_CM090947.1.recode.vcf
rm noTexas_CM090672.1.recode.vcf

rm desert_CM090947.1.recode.vcf
rm desert_CM090672.1.recode.vcf

# cat together
grep "^#" noTexas_CM090671.1.recode.vcf > noTexas.vcf

for i in $( ls noTexas_*.recode.vcf ); do grep -v "^#" $i >> noTexas.vcf; done

grep "^#" desert_CM090671.1.recode.vcf > desert.vcf

for i in $( ls desert_*.recode.vcf ); do grep -v "^#" $i >> desert.vcf; done

# remove individual scaffolds
rm *recode*

########################################
######## PCA 1 and ADMIXTURE
########################################

# pca for desert bighorn no Texas
# make chromosome map for the vcf
grep -v "#" noTexas.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the vcf intended for PCA with plink output
vcftools --vcf noTexas.vcf --plink --chrom-map chrom_map.txt --out pca_noTexas

# convert with plink for pca dataset
plink --file pca_noTexas --recode12 --allow-extra-chr --out pca_noTexas_plink

# run  dataset for pca
plink --file pca_noTexas_plink --pca --allow-extra-chr --out pca_noTexas_plink_pca


# admixture
# run vcftools for the vcf intended for admixture with plink output
vcftools --vcf noTexas.vcf --plink --chrom-map chrom_map.txt --out admixture_noTexas

# convert with plink for admixture dataset
plink --file admixture_noTexas --recode12 --allow-extra-chr --out admixture_noTexas_plink

# run admixture
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv admixture_noTexas_plink.ped $K  | tee log_${K}.out; done

# look at logs
cat log_*out | grep "CV error" # lowest CV error = (K=4): 0.42774

########################################
######## PCA 2 and ADMIXTURE
########################################

# pca for all desert bighorn
# make chromosome map for the vcf
grep -v "#" desert.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the vcf intended for PCA with plink output
vcftools --vcf desert.vcf --plink --chrom-map chrom_map.txt --out pca_desert

# convert with plink for pca dataset
plink --file pca_desert --recode12 --allow-extra-chr --out pca_desert_plink

# run  dataset for pca
plink --file pca_desert_plink --pca --allow-extra-chr --out pca_desert_plink_pca


# admixture
# run vcftools for the vcf intended for admixture with plink output
vcftools --vcf desert.vcf --plink --chrom-map chrom_map.txt --out admixture_desert

# convert with plink for admixture dataset
plink --file admixture_desert --recode12 --allow-extra-chr --out admixture_desert_plink


# run admixture
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv admixture_desert_plink.ped $K  | tee log_${K}.out; done

# look at logs
cat log_*out | grep "CV error" # lowest CV error = (K=6): 0.32584





