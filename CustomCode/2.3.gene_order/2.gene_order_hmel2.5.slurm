#!/bin/bash
#SBATCH -n 12 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 600 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=24000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o gene_order_%j.out # Standard out goes to this file
#SBATCH -e gene_order_%j.err # Standard err goes to this filehostname



## -------------------------------------------------------------------------------
## map genes to MEDUSA assembly :: NUCMER :: -------------------------------------

# parameters
species="/n/scratchlfs/mallet_lab/fseixas/medusa/w2.100gap/hele2hmel/hele2hmel_min05kb.100gap.w2.fasta"
genes="/n/scratchlfs/mallet_lab/fseixas/gene_order/to_hmel/hmel2.5.genes.buf500.bp.fa"
output="hele2hmel.m05kb.w2.fixgap.hmelgenes500bp"
newfas=`echo $species | sed 's/.fasta/.fa/' | awk -F'[/]' '{print $NF}' `

# dirs
mkdir hele.m05kb.w2.fixgap
cd hele.m05kb.w2.fixgap

# load modules
module load mummer/4.0.0beta2-fasrc01

# create symbolic links
ln -s $species $newfas

# map genes to MEDUSA assembly
nucmer --threads=12 --prefix=$output $newfas $genes
show-coords -Bd $output.delta > $output.coords



## ----------------------------------------------------------------
## get best mappings (in R) ---------------------------------------
module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.1.nucmer2tophits.R -i $output.coords -o $output -b 2500 -c 12


## ----------------------------------------------------------------
## make summary of gene hits per chromosome -----------------------
# for each hmel chromosome count hits of respective genes in MEDUSA scaffolds 
Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.2.scaffolds_per_chromosome.R -g ../hmel2.5.genes_coord.tab -t hele2hmel.m05kb.w2.fixgap.hmelgenes500bp.tophits.txt -p hele2hmel.m05kb.w2.fixgap.hmelgenes500bp

# test correlation of gene coordinates between reference (e.g. Hmel2.5) and MEDUSA assemblies [testing if e.g. genes from chr01 mapping to the respective scaffold in the MEDUSA assembly are in the expected order] 
Rscript /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/2.3.gene_coord_correlations.R -g ../hmel2.5.genes_coord.tab -t $output.tophits.txt -p $output



