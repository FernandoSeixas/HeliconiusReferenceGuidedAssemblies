#!/bin/bash
#SBATCH -n 8    # Number of cores requested
#SBATCH -N 1    # Ensure that all cores are on one machine
#SBATCH -t 120 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=48000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o logs/liftoff_%j.out # Standard out goes to this file
#SBATCH -e logs/liftoff_%j.err # Standard err goes to this filehostname


## variables
ref=$1
spp=$2


## load modules
module load Anaconda3/5.1.0
source activate liftoff


## create directory to paste reference for this particularly run 
mkdir run.$spp.$ref; 

if [ $ref = "hmelv25" ]; then 
    cp 0.references/Hmel2.5.gff3 run.$spp.$ref/$ref.gff3;
    cp 0.references/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa run.$spp.$ref/$ref.scaffolds.fa;
    cp 0.finalAssemblies/$spp-2-hmelv25.HM.100gap.finalAssembly.fasta run.$spp.$ref;
fi

if [ $ref = "heradem" ]; then 
    cp 0.references/Heliconius_erato_v1_-_genes.final.gff run.$spp.$ref/$ref.gff3;
    #cp 0.references/unused/Heliconius_erato_demophoon_v1.gff3 run.$spp.$ref/$ref.gff3;
    cp 0.references/Heliconius_erato_demophoon_v1_-_scaffolds.fa run.$spp.$ref/$ref.scaffolds.fa;
    cp 0.finalAssemblies/$spp-2-heradem.HM.100gap.finalAssembly.fasta run.$spp.$ref;
fi


## run liftoff
liftoff \
    -p 8 \
    -g run.$spp.$ref/$ref.gff3 \
    -chroms 0.scaffAssociations/$spp-2-$ref.bestAssociation.txt \
    -unplaced 0.scaffAssociations/$spp-2-$ref.Unplaced.txt \
    -d 2 \
    -flank 0.2 \
    -o $spp-2-$ref.HM.100gap.finalAssembly.gff \
    -u $spp-2-$ref.unmapped.txt \
    -dir temp.$spp.$ref \
    run.$spp.$ref/$spp-2-$ref.HM.100gap.finalAssembly.fasta \
    run.$spp.$ref/$ref.scaffolds.fa


