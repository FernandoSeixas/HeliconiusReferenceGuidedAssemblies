# load modules
module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER

#variables
ref=$1

# symbolic link to original genome
if [ "$ref" == "heradem" ]; then ln -s /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hera.dem.v1/Heliconius_erato_demophoon_v1_-_scaffolds.fa heradem.fa; fi
if [ "$ref" == "hmelv25" ]; then ln -s /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmel2.5/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa hmelv25.fa; fi

# get scaffold sizes
awk '{ if ($1 ~ /^>/) printf "%s\t", $1; else print length($1)}' $ref.fa | sed 's,>,,g' | grep -v "Herato_mt" > $ref.scaffold.sizes

# scaffold 2 chromosome
perl scaffs2chroms.pl $ref.scaffold.sizes $ref > $ref.scaffGroups

# update gene coordinates
Rscript genecoord.vGlued.R
