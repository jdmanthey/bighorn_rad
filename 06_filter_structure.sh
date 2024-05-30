# interactive session with 4 cpu and 32 GB memory

source activate bcftools

workdir=/lustre/scratch/jmanthey/10_bighorn_rad/

cd ${workdir}
cd 22_filter_script

########################################
######## filter for first 8 datasets
########################################

# filter for dataset 1 
# for RAxML
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset1.txt \
--max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset1

# filter for dataset 2
# for observed heterozygosity, genetic distance, splitstree
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset2.txt \
--max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset2

# filter for dataset 3
# for PCA outside Texas
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset3.txt \
--min-alleles 2 --max-alleles 2 --mac 3 --thin 10000 --max-maf 0.49 \
--max-missing 0.94 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset3

# filter for dataset 4
# for PCA, ADMIXTURE all desert bighorn
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset4.txt \
--min-alleles 2 --max-alleles 2 --mac 3 --thin 10000 --max-maf 0.49 \
--max-missing 0.953 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset4

# filter for dataset 5
# for kinship, EEMS all Texas samples
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset5.txt \
--min-alleles 2 --max-alleles 2 --mac 3 --thin 10000 --max-maf 0.49 \
--max-missing 1.0 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset5

# filter for dataset 6
# for EEMS region 1
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset6.txt \
--min-alleles 2 --max-alleles 2 --mac 3 --thin 10000 --max-maf 0.49 \
--max-missing 1.0 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset6

# filter for dataset 7
# for EEMS region 2
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset7.txt \
--min-alleles 2 --max-alleles 2 --mac 3 --thin 10000 --max-maf 0.49 \
--max-missing 1.0 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset7

# filter for dataset 8
# for pairwise Fst of all populations with at least 4 individuals
vcftools --vcf ${workdir}/05_filtered_vcf/total2.recode.vcf --keep keep_dataset8.txt \
--min-alleles 2 --max-alleles 2 --max-maf 0.49 \
--recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/dataset8


########################################
######## Copy filtered VCFs to analysis directories
########################################

cd ../05_filtered_vcf

cp dataset1.recode.vcf ../07_phylogeny
cp dataset2.recode.vcf ../07_phylogeny
cp dataset2.recode.vcf ../08_stats
cp dataset3.recode.vcf ../06_structure
cp dataset4.recode.vcf ../06_structure
cp dataset5.recode.vcf ../09_relatedness
cp dataset5.recode.vcf ../11_eems
cp dataset6.recode.vcf ../11_eems
cp dataset7.recode.vcf ../11_eems
cp dataset8.recode.vcf ../08_stats

########################################
######## PCA 1
########################################

cd ../06_structure

# pca for desert bighorn no Texas
# make chromosome map for the vcf
grep -v "#" dataset3.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the vcf intended for PCA with plink output
vcftools --vcf dataset3.recode.vcf --plink --chrom-map chrom_map.txt --out pca_dataset3

# convert with plink for pca dataset
plink --file pca_dataset3 --recode12 --allow-extra-chr --out pca_dataset3_plink

# run  dataset for pca
plink --file pca_dataset3_plink --pca --allow-extra-chr --out pca_dataset3_plink_pca

########################################
######## PCA 2 and ADMIXTURE
########################################

# pca for desert bighorn
# make chromosome map for the vcf
grep -v "#" dataset4.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the vcf intended for PCA with plink output
vcftools --vcf dataset4.recode.vcf --plink --chrom-map chrom_map.txt --out pca_dataset4

# convert with plink for pca dataset
plink --file pca_dataset4 --recode12 --allow-extra-chr --out pca_dataset4_plink

# run  dataset for pca
plink --file pca_dataset4_plink --pca --allow-extra-chr --out pca_dataset4_plink_pca


# admixture
# run vcftools for the vcf intended for admixture with plink output
vcftools --vcf dataset4.recode.vcf --plink --chrom-map chrom_map.txt --out admixture

# convert with plink for admixture dataset
plink --file admixture --recode12 --allow-extra-chr --out admixture_plink


# run admixture
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv admixture_plink.ped $K  | tee log_${K}.out; done

########################################
######## Fst
########################################

cd ../08_stats
cp ../22_filter_script/keep_dataset8.txt .

# get population 'keep' files for vcftools fst
grep "^Arizona_" keep_dataset8.txt > Arizona_fst_.txt
grep "^GabbsValley-GBL_" keep_dataset8.txt > GabbsValley-GBL_fst_.txt
grep "^GabbsValley-ML_" keep_dataset8.txt > GabbsValley-ML_fst_.txt
grep "^MuddyMtnsNV_" keep_dataset8.txt > MuddyMtnsNV_fst_.txt
grep "^MormonMtnsNV_" keep_dataset8.txt > MormonMtnsNV_fst_.txt
grep "^BlackMtnsAZ_" keep_dataset8.txt > BlackMtnsAZ_fst_.txt
grep "^SonoraMX_" keep_dataset8.txt > SonoraMX_fst_.txt
grep "^SierraDiablo_" keep_dataset8.txt > SierraDiablo_fst_.txt
grep "^Beaches_" keep_dataset8.txt > Beaches_fst_.txt
grep "^Baylor_" keep_dataset8.txt > Baylor_fst_.txt
grep "^VanHorn_" keep_dataset8.txt > VanHorn_fst_.txt
grep "^ElephantMtn_" keep_dataset8.txt > ElephantMtn_fst_.txt
grep "^Dove_" keep_dataset8.txt > Dove_fst_.txt
grep "^Basse_" keep_dataset8.txt > Basse_fst_.txt
grep "^BlackGap_" keep_dataset8.txt > BlackGap_fst_.txt
grep "^SouthernCalifornia_" keep_dataset8.txt > SouthernCalifornia_fst_.txt
grep "^RockyMtnNM_" keep_dataset8.txt > RockyMtnNM_fst_.txt



# nested loop for pairwise comparisons to get fst
echo "pop1 pop2 fst" > fst_output.txt
for i in *_fst_.txt
do
  for j in *_fst_.txt
  do
    if [ "$i" \< "$j" ]
    then
     vcftools --vcf ${workdir}/08_stats/dataset8.recode.vcf --weir-fst-pop $i --weir-fst-pop $j \
     --out ${i%_fst_.txt}__${j%_fst_.txt}
     fst=$( grep "mean" ${i%_fst_.txt}__${j%_fst_.txt}.log | cut -d ' ' -f 7 )
     pop1=${i%_fst_.txt}
     pop2=${j%_fst_.txt}
     echo "$pop1 $pop2 $fst" >> fst_output.txt
    fi
  done
done

# remove unneeded files
rm *_fst_.txt
rm *.weir.fst
rm *log





