# move to working directory


# index the reference
interactive -p nocona

cd references

source activate bcftools

bwa-mem2 index GCA_042477335.2_ARS-UI_OviCan_v2_genomic.fna 

samtools faidx GCA_042477335.2_ARS-UI_OviCan_v2_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCA_042477335.2_ARS-UI_OviCan_v2_genomic.fna \
O=/home/jmanthey/references/GCA_042477335.2_ARS-UI_OviCan_v2_genomic.dict





# make directories for organization during analyses
mkdir 00_fastq
mkdir 01_cleaned
mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 04_vcf
mkdir 05_filtered_vcf
mkdir 06_structure
mkdir 07_phylogeny
mkdir 08_stats
mkdir 09_relatedness
mkdir 10_popsize
mkdir 11_eems
mkdir 20_align_script
mkdir 21_filter_script

# put all fastq files in the 00_fastq directory
# rename all samples
cd 00_fastq 

while read -r name1 name2; do
	mv $name1 $name2
done < rename.txt
