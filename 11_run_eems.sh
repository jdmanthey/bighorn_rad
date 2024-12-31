#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=eems
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate eems

/home/jmanthey/eems/runeems_snps/src/runeems_snps --params Texas_bighorn_eems.params --seed 123
