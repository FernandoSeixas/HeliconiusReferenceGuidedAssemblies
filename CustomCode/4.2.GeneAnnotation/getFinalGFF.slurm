#!/bin/bash
#SBATCH -n 1    # Number of cores requested
#SBATCH -N 1    # Ensure that all cores are on one machine
#SBATCH -t 60 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=8000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o logs/finalgff_%j.out # Standard out goes to this file
#SBATCH -e logs/finalgff_%j.err # Standard err goes to this filehostname

## variables
ref=$1
spp=$2
export ref=$ref
export spp=$spp

## load modules
module load perl/5.10.1-fasrc05
module load perl-modules/5.10.1-fasrc13

## extract valid transcripts and corresponding annotations
rm $spp-2-$ref.tmp.gff
touch $spp-2-$ref.tmp.gff
# go through each gene
for gn in `cut -f1 $spp-2-$ref.TranscriptsValidORFs.txt | sort | uniq `; do
    # write gene feature to file
    grep -w $gn 1.gff/v3/$spp-2-$ref.HM.100gap.finalAssembly.gff | grep "gene" >> $spp-2-$ref.tmp.gff
    # go through each transcript and write mRNA, exon and CDS features to file
    for id in `grep -w $gn $spp-2-$ref.TranscriptsValidORFs.txt | cut -f2 `; do
        tid=`echo "$id;extra"`
        grep -Fr $tid 1.gff/v3/$spp-2-$ref.HM.100gap.finalAssembly.gff >> $spp-2-$ref.tmp.gff
    done
done

## sort gff file
module load Anaconda/5.0.1-fasrc02
source activate agat

agat_convert_sp_gxf2gxf.pl --gff $spp-2-$ref.tmp.gff --output $spp-2-$ref.final.gff


