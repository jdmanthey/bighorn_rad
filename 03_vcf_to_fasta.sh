interactive -p nocona -c 4

source activate bcftools

workdir=/lustre/scratch/jmanthey/07_bighorn_rad

cd ${workdir}/07_phylogeny

# remove sex chromosomes
rm phylo_CM090947.1.recode.vcf
rm phylo_CM090672.1.recode.vcf

##########################################
### cat all individual scaffolds together
##########################################

grep "^#" phylo_CM090671.1.recode.vcf > phylo.vcf

for i in $( ls phylo_*.recode.vcf ); do grep -v "^#" $i >> phylo.vcf; done

rm *recode*

##########################################
### convert to FASTA
##########################################

# do this in R
options(scipen=999)

# read in a small vcf (don't use for large vcf files)
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}

# read in vcf  
x <- read_vcf("phylo.vcf")

output_names_fasta <- paste(">", colnames(x[,10:ncol(x)]), sep="")

# subset the genotypes from the allele info
allele_info <- x[,4:5]
genotypes <- x[,10:ncol(x)]

# keep only first three characters for each genotype column
for(a in 1:ncol(genotypes)) {
	genotypes[,a] <- substr(genotypes[,a], 1, 3)
}

# convert all numbers in genotypes to actual bases, heterozygous sites as random allele draw of the two
for(a in 1:nrow(genotypes)) {
	if(allele_info[a,2] == ".") { # if non-polymorphic
        genotypes[a,] <- gsub("0/0", allele_info[a,1], genotypes[a,])
        genotypes[a,] <- gsub("\\./\\.", "?", genotypes[a,])
    } else { # if polymorphic
        both_genotypes <- sort(as.character(allele_info[a,]))
        if(both_genotypes[1] == "A" & both_genotypes[2] == "C") { het = c("A", "C") }
        if(both_genotypes[1] == "A" & both_genotypes[2] == "G") { het = c("A", "G") }
        if(both_genotypes[1] == "A" & both_genotypes[2] == "T") { het = c("A", "T") }
        if(both_genotypes[1] == "C" & both_genotypes[2] == "G") { het = c("C", "G") }
        if(both_genotypes[1] == "C" & both_genotypes[2] == "T") { het = c("C", "T") }
        if(both_genotypes[1] == "G" & both_genotypes[2] == "T") { het = c("G", "T") }
        genotypes[a,] <- gsub("0/0", allele_info[a,1], genotypes[a,])
        genotypes[a,] <- gsub("\\./\\.", "?", genotypes[a,])
        genotypes[a,] <- gsub("1/1", allele_info[a,2], genotypes[a,])
        genotypes[a,][genotypes[a,] == "0/1"] <- sample(het, length(genotypes[a,][genotypes[a,] == "0/1"]), replace=T)
    }
}

# remove any sites that are invariant (due to random sampling of 
# variant sites only found in heterozygous state)
keep <- list()
for(a in 1:nrow(genotypes)) {
	a_rep <- genotypes[a,]
	a_rep <- a_rep[a_rep != "?"]
	if(length(unique(a_rep)) > 1) {
		keep[[a]] <- a
	}
}
keep <- unlist(keep)
genotypes <- genotypes[keep,]


# output name
outname <- "bighorn.fasta"

# write output
for(a in 1:ncol(genotypes)) {
    if(a == 1) {
    	write(output_names_fasta[a], file=outname, ncolumns=1)
        write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
    } else {
        write(output_names_fasta[a], file=outname, ncolumns=1, append=T)
        write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
	}
}
dim(genotypes)
