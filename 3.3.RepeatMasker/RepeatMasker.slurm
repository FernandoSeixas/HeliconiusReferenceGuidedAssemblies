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
spack load perl
module load perl-5.36.0-gcc-12.2.0-2qh2sr4
module load perl_modules/5.10.1-ncf

## variables
spp=$1
export spp=$spp


## filter only large scaffolds in the final assemblies 
awk -v min=1000000 'length($1) > min {printf "%s\n%s\n",f,$1} {f=$1}' \
    /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/finalAssemblies/hmelv25/$spp-2-hmelv25.HM.100gap.finalAssembly.fasta \
    > $spp-2-hmelv25.HM.100gap.bigScaffs.fasta


## run RepeatMasker
~/software/RepeatMasker/RepeatMasker \
    -pa 24 \
    -species Heliconius \
    -dir /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/$spp.bigScaffs.Heliconius \
    -s \
    -xsmall \
    -nocut \
    $spp-2-hmelv25.HM.100gap.bigScaffs.fasta


## Parser Output to generate landscape of Repeats through time [%div]
perl /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/Parsing-RepeatMasker-Outputs/parseRM.pl \
    -i /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/$spp.bigScaffs.Heliconius/$spp-2-hmelv25.HM.100gap.bigScaffs.fasta.out -p -l 50,1
