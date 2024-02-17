#!/bin/bash
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 10 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o gene_order_%j.out # Standard out goes to this file
#SBATCH -e gene_order_%j.err # Standard err goes to this filehostname

# parameters
input=$1
genes=$2
refnam=$3
bfsize=$4
minlen=$5
prefix=${input##*/}
spname=${prefix/.fasta/}
spp=${spname/\-2\-hmelv25.HM.100gap.finalAssembly/}
output=$spname.$refnam.genes.$bfsize-bp

# create directories
mkdir $prefix
cd $prefix
echo $spname

### -------------------------------------------------------------------------------
#module load mummer/4.0.0beta2-fasrc01
## map genes to MEDUSA assembly :: NUCMER :: -------------------------------------
#nucmer --threads=8 --prefix=$output $input ../$genes
#show-coords -Bd $output.delta > $output.coords


## ----------------------------------------------------------------
module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
## get best mappings (in R) ---------------------------------------
#Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.1.nucmer2tophits.R -i $output.coords -o $output -b $bfsize -c 8
## ----------------------------------------------------------------
## make summary of gene hits per chromosome -----------------------
# for each hmel chromosome count hits of respective genes in MEDUSA scaffolds 
#Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.2.scaffolds_per_chromosome.R -g ~/code/heliconius_seixas/1.pseudo_references/0.1.references/$refnam.genes_coord.vGlued.tab -t $output.tophits.txt -p $output.mgl$minlen -m $minlen
#Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.2.scaffolds_per_chromosome.v2.R -s $spp -r $refnam -g ~/code/heliconius_seixas/1.pseudo_references/0.1.references/$refnam.genes_coord.vOrigin.tab -t $output.tophits.txt -p $output.mgl$minlen -m $minlen
## test correlation of gene coordinates between reference (e.g. Hmel2.5) and MEDUSA assemblies [testing if e.g. genes from chr01 mapping to the respective scaffold in the MEDUSA assembly are in the expected order] 
Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.3.gene_coord_correlations.v2.R \
-r $refnam \
-s $spp \
-g ~/code/heliconius_seixas/1.pseudo_references/0.1.references/$refnam.genes_coord.vOrigin.tab \
-t $output.tophits.txt \
-p $output.mgl$minlen \
-m $minlen


