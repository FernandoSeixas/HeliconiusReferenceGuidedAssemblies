# load modules
module load minimap2/2.9-fasrc01


## variables ===================
# Original References
#REF="hmelv25"
#TARGET="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa"
#TARGET="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/hmelv25/hmelv25.glue.fa"
REF="heradem"
TARGET="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/Heliconius_erato_demophoon_v1_-_scaffolds.fa"
#TARGET="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/heradem/heradem.glue.fa"
#REF="heralat"
#TARGET="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hera.lat.v1/Heliconius_erato_lativitta_v1_-_scaffolds.fa"
#TARGET="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/heralat/heralat.glue.fa"
# Glued
#REF="hmelv25Glued
#TARGET="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/hmelv25/hmelv25.glue.fa"
#REF="herademGlued"
#TARGET="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/heradem/heradem.glue.fa"
#REF="heralatGlued"
#TARGET="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/heralat/heralat.glue.fa"
# 
MINSIZE=5000

## create directories ==================================
mkdir 1.subsetContigs
mkdir 2.minimap2


## map DISCOVAR/w2rap to references ===================
for SPP in `cat ~/code/heliconius_seixas/1.pseudo_references/spp.list`; do
    QUERY="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/2.1.medusa/vReferGenomes/assemblies/"$SPP"_DISCOVAR_cleaned_A_ref.sline.fa"
    export SPP=$SPP
    export REF=$REF
    export QUERY=$QUERY
    export TARGET=$TARGET
    # filter big DISCOVAR contigs
    awk -v min=5000 'length($1) >= min {printf "%s\n%s\n",f,$1} {f=$1}' $QUERY > 1.subsetContigs/$SPP.scontigs-5kb.fasta
    # map scaffolds to assembly
    minimap2 $TARGET 1.subsetContigs/$SPP.scontigs-5kb.fasta > 2.minimap2/$SPP-2-$REF.mm2
    # get primary alignments and filter based on length
    grep "tp:A:P" 2.minimap2/$SPP-2-$REF.mm2 | awk '{if ($11 > 1000) print }' > 2.minimap2/$SPP-2-$REF.filter.mm2
done


module load minimap2/2.9-fasrc01

#REF="hmelv25"
#TARGET="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa"
REF="heradem"
TARGET="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/Heliconius_erato_demophoon_v1_-_scaffolds.fa"
MINSIZE=5000

SPP="etal_w2rap"
QUERY="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/DISCOCAR_w2rap/etal.DISCOVAR.scaffolds.fa"

export SPP=$SPP
export REF=$REF
export QUERY=$QUERY
export TARGET=$TARGET
# filter big DISCOVAR contigs
awk -v min=5000 'length($1) >= min {printf "%s\n%s\n",f,$1} {f=$1}' $QUERY > 1.subsetContigs/$SPP.scontigs-5kb.fasta
# map scaffolds to assembly
minimap2 $TARGET 1.subsetContigs/$SPP.scontigs-5kb.fasta > 2.minimap2/$SPP-2-$REF.mm2
# get primary alignments and filter based on length
grep "tp:A:P" 2.minimap2/$SPP-2-$REF.mm2 | awk '{if ($11 > 1000) print }' > 2.minimap2/$SPP-2-$REF.filter.mm2
