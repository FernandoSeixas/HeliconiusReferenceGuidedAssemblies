#!/bin/bash
#SBATCH -n 24       # Number of cores requested
#SBATCH -N 1        # Ensure that all cores are on one machine
#SBATCH -t 7200     # Runtime in minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o RepeatMasker_%j.out # Standard out goes to this file
#SBATCH -e RepeatMasker_%j.err # Standard err goes to this filehostname


## load modules
module load perl/5.26.1-fasrc01
module load perl-modules/5.26.1-fasrc01


## filter only large scaffolds in the final assemblies
awk -v min=1000000 'length($1) > min {printf "%s\n%s\n",f,$1} {f=$1}' \
    /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmel2.5/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa \
    > hmelv25.bigScaffs.fasta


## run RepeatMasker
~/software/RepeatMasker/RepeatMasker \
    -pa 24 \
    -species Heliconius \
    -dir /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/hmelv25.bigScaffs.Heliconius \
    -a \
    -s \
    -xsmall \
    -nocut \
    hmelv25.bigScaffs.fasta


## Parser Output to generate landscape of Repeats through time [%div]
perl /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/Parsing-RepeatMasker-Outputs/parseRM.pl \
    -i /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out -p -l 50,1
