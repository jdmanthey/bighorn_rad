options(scipen=999)

x_stats <- list.files("output_files", pattern="*stats.txt", full.names=T)

# summarize stats
fis <- list()
het <- list()
pi <- list()
Dxy <- list()
Fst <- list()
private <- list()
shared <- list()
fixed <- list()
for(a in 1:length(x_stats)) {
	a_rep <- read.table(x_stats[a], sep="\t", header=T)
	fis[[a]] <- a_rep[a_rep[,3] == "fis",]
	het[[a]] <- a_rep[a_rep[,3] == "heterozygosity",]
	pi[[a]] <- a_rep[a_rep[,3] == "pi",]
	Dxy[[a]] <- a_rep[a_rep[,3] == "Dxy",]
	Fst[[a]] <- a_rep[a_rep[,3] == "Fst",]
	private[[a]] <- a_rep[a_rep[,3] == "private",]
	shared[[a]] <- a_rep[a_rep[,3] == "shared",]
	fixed[[a]] <- a_rep[a_rep[,3] == "fixed",]
}

fis <- do.call(rbind, fis)
het <- do.call(rbind, het)
pi <- do.call(rbind, pi)
Dxy <- do.call(rbind, Dxy)
Fst <- do.call(rbind, Fst)
private <- do.call(rbind, private)
shared <- do.call(rbind, shared)
fixed <- do.call(rbind, fixed)

# Fis
output_rep <- c()
a_stat_matrix <- fis
for(a in 1:length(unique(a_stat_matrix[,1]))) {
	a_ID <- unique(a_stat_matrix[,1])[a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_ID,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat * a_rep$number_sites) / a_num_sites
	output_rep <- rbind(output_rep, c(a_ID, a_num_sites, a_num_var_sites, a_stat))
}
fis_output <- data.frame(ID=as.character(output_rep[,1]), number_sites=as.numeric(output_rep[,2]), number_variable_sites=as.numeric(output_rep[,3]), calculated_statistic=as.numeric(output_rep[,4]))


# heterozygosity
output_rep <- c()
a_stat_matrix <- het
for(a in 1:length(unique(a_stat_matrix[,1]))) {
	a_ID <- unique(a_stat_matrix[,1])[a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_ID,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat * a_rep$number_sites) / a_num_sites
	output_rep <- rbind(output_rep, c(a_ID, a_num_sites, a_num_var_sites, a_stat))
}
het_output <- data.frame(ID=as.character(output_rep[,1]), number_sites=as.numeric(output_rep[,2]), number_variable_sites=as.numeric(output_rep[,3]), calculated_statistic=as.numeric(output_rep[,4]))


# pi
output_rep <- c()
a_stat_matrix <- pi
for(a in 1:length(unique(a_stat_matrix[,1]))) {
	a_ID <- unique(a_stat_matrix[,1])[a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_ID,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat * a_rep$number_sites) / a_num_sites
	output_rep <- rbind(output_rep, c(a_ID, a_num_sites, a_num_var_sites, a_stat))
}
pi_output <- data.frame(ID=as.character(output_rep[,1]), number_sites=as.numeric(output_rep[,2]), number_variable_sites=as.numeric(output_rep[,3]), calculated_statistic=as.numeric(output_rep[,4]))


# private
output_rep <- c()
a_stat_matrix <- private
for(a in 1:length(unique(a_stat_matrix[,1]))) {
	a_ID <- unique(a_stat_matrix[,1])[a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_ID,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat)
	output_rep <- rbind(output_rep, c(a_ID, a_num_sites, a_num_var_sites, a_stat))
}
private_output <- data.frame(ID=as.character(output_rep[,1]), number_sites=as.numeric(output_rep[,2]), number_variable_sites=as.numeric(output_rep[,3]), calculated_statistic=as.numeric(output_rep[,4]))



# shared
output_rep <- c()
a_stat_matrix <- shared
for(a in 1:length(unique(a_stat_matrix[,1]))) {
	a_ID <- unique(a_stat_matrix[,1])[a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_ID,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat)
	output_rep <- rbind(output_rep, c(a_ID, a_num_sites, a_num_var_sites, a_stat))
}
shared_output <- data.frame(ID=as.character(output_rep[,1]), number_sites=as.numeric(output_rep[,2]), number_variable_sites=as.numeric(output_rep[,3]), calculated_statistic=as.numeric(output_rep[,4]))


# fixed
output_rep <- c()
a_stat_matrix <- fixed
for(a in 1:length(unique(a_stat_matrix[,1]))) {
	a_ID <- unique(a_stat_matrix[,1])[a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_ID,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat)
	output_rep <- rbind(output_rep, c(a_ID, a_num_sites, a_num_var_sites, a_stat))
}
fixed_output <- data.frame(ID=as.character(output_rep[,1]), number_sites=as.numeric(output_rep[,2]), number_variable_sites=as.numeric(output_rep[,3]), calculated_statistic=as.numeric(output_rep[,4]))


# fst
output_rep <- c()
a_stat_matrix <- Fst
populations <- unique(c(a_stat_matrix[,1], a_stat_matrix[,2]))
pop_combs <- combn(populations,2)
for(a in 1:ncol(pop_combs)) {
	a_pop1 <- pop_combs[1,a]
	a_pop2 <- pop_combs[2,a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_pop1 & a_stat_matrix[,2] == a_pop2,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat * a_rep$number_variable_sites) / a_num_var_sites
	output_rep <- rbind(output_rep, c(a_pop1, a_pop2, a_num_sites, a_num_var_sites, a_stat))
}
fst_output <- data.frame(Pop1=as.character(output_rep[,1]), Pop2=as.character(output_rep[,2]), number_sites=as.numeric(output_rep[,3]), number_variable_sites=as.numeric(output_rep[,4]), calculated_statistic=as.numeric(output_rep[,5]))


# dxy
output_rep <- c()
a_stat_matrix <- Dxy
populations <- unique(c(a_stat_matrix[,1], a_stat_matrix[,2]))
pop_combs <- combn(populations,2)
for(a in 1:ncol(pop_combs)) {
	a_pop1 <- pop_combs[1,a]
	a_pop2 <- pop_combs[2,a]
	a_rep <- a_stat_matrix[a_stat_matrix[,1] == a_pop1 & a_stat_matrix[,2] == a_pop2,]
	a_num_sites <- sum(a_rep$number_sites)
	a_num_var_sites <- sum(a_rep$number_variable_sites)
	a_stat <- sum(a_rep$calculated_stat * a_rep$number_variable_sites) / a_num_var_sites
	output_rep <- rbind(output_rep, c(a_pop1, a_pop2, a_num_sites, a_num_var_sites, a_stat))
}
dxy_output <- data.frame(Pop1=as.character(output_rep[,1]), Pop2=as.character(output_rep[,2]), number_sites=as.numeric(output_rep[,3]), number_variable_sites=as.numeric(output_rep[,4]), calculated_statistic=as.numeric(output_rep[,5]))


# write output tables
write.table(fis_output, file="_stats_fis.txt", sep="\t", quote=F, row.names=F)
write.table(het_output, file="_stats_het.txt", sep="\t", quote=F, row.names=F)
write.table(pi_output, file="_stats_pi.txt", sep="\t", quote=F, row.names=F)
write.table(shared_output, file="_stats_shared.txt", sep="\t", quote=F, row.names=F)
write.table(private_output, file="_stats_private.txt", sep="\t", quote=F, row.names=F)
write.table(fixed_output, file="_stats_fixed.txt", sep="\t", quote=F, row.names=F)
write.table(fst_output, file="_stats_fst.txt", sep="\t", quote=F, row.names=F)
write.table(dxy_output, file="_stats_dxy.txt", sep="\t", quote=F, row.names=F)






# summarize distance matrices
x_diffs <- list.files("output_files", pattern="*diffs", full.names=T)

number_sites <- c()
output <- c()
for(a in 1:length(x_diffs)) {
	a_sites <- as.numeric(strsplit(x_diffs[a], "__")[[1]][[2]])
	a_rep <- read.table(x_diffs[a], header=T)
	a_rep <- a_rep * a_sites
	if(a == 1) {
		output <- a_rep
	} else {
		output <- output + a_rep
	}
	number_sites <- c(number_sites, a_sites)
}
number_sites <- sum(number_sites)
output_standardized <- output / number_sites
rownames(output_standardized) <- colnames(output_standardized)

write.table(output_standardized, file="_stat_distance_matrix.txt", sep="\t", quote=F)


# check okay format of output matrix
test <- read.table("_stat_distance_matrix.txt", header=T, row.names=1)



