#!/bin/bash
#SBATCH -J haplomerger  # A single job name for the array
#SBATCH -n 8   # Number of cores requested
#SBATCH -N 1    # Ensure that all cores are on one machine
#SBATCH -t 1200 # Runtime in minutes
#SBATCH -p shared  # Partition to submit to
#SBATCH --mem=32000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o haploM_%A.out # Standard out goes to this file
#SBATCH -e haploM_%A.err # Standard err goes to this filehostname

# load modules and add dirs to PATH
module load perl/5.26.1-fasrc01
module load lastz/1.04.00-fasrc01
module load ucsc/20150820-fasrc01
PATH=/n/home12/fseixas/software/assemblies/HaploMerger2_20180603/chainNet_jksrc20100603_centOS6/:/n/home12/fseixas/software/assemblies/HaploMerger2_20180603/gapCloser_v1.12/:/n/home12/fseixas/software/assemblies/HaploMerger2_20180603/SSPACE-STANDARD-3.0_linux-x86_64/:$PATH
echo $PATH

# variables
SPP=$1

# run HaploMerger
bash run_all.batch

# convert assemblies to single line fasta
zcat $SPP\_DISCOVAR_cleaned.fa.gz | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > $SPP\_DISCOVAR_cleaned.sline.fa
bash ~/code/heliconius_seixas/1.pseudo_references/2.1.medusa/assembly_summarySimple.sh -i $SPP\_DISCOVAR_cleaned_ref.sline.fa -n $SPP
zcat $SPP\_DISCOVAR_cleaned_A.fa.gz | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > $SPP\_DISCOVAR_cleaned_A.sline.fa
bash ~/code/heliconius_seixas/1.pseudo_references/2.1.medusa/assembly_summarySimple.sh -i $SPP\_DISCOVAR_cleaned_A.sline.fa -n $SPP
zcat $SPP\_DISCOVAR_cleaned_A_ref.fa.gz | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > $SPP\_DISCOVAR_cleaned_A_ref.sline.fa
bash ~/code/heliconius_seixas/1.pseudo_references/2.1.medusa/assembly_summarySimple.sh -i $SPP\_DISCOVAR_cleaned_A_ref.sline.fa -n $SPP

