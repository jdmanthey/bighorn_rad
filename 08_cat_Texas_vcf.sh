# interactive session with 4 cpu and 32 GB memory

interactive -p nocona -c 4 -m 8G

source activate bcftools

workdir=/lustre/scratch/jmanthey/07_bighorn_rad/

cd ${workdir}/09_relatedness

########################################
######## combine files
########################################

# remove sex chromosomes
rm Texas_CM090947.1.recode.vcf
rm Texas_CM090672.1.recode.vcf

# cat together
grep "^#" Texas_CM090671.1.recode.vcf > Texas.vcf

for i in $( ls Texas_*.recode.vcf ); do grep -v "^#" $i >> Texas.vcf; done

# remove individual scaffolds
rm *recode*


