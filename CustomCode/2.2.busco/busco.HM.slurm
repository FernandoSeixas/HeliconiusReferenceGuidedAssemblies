#!/bin/bash
#SBATCH -n 8   # Number of cores requested
#SBATCH -N 1    # Ensure that all cores are on one machine
#SBATCH -t 600 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=16000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o busco_%j.out # Standard out goes to this file
#SBATCH -e busco_%j.err # Standard err goes to this filehostname

## define variables
SPP=$1
FAS="${SPP}_DISCOVAR_cleaned_A_ref.sline.fa"


# state what you're doing to console
echo "running busco for $SPP"

## load modules ========================================
module load blast/2.6.0+-fasrc01
module load augustus/3.3-fasrc02
module load hmmer/3.1b2-fasrc01
export AUGUSTUS_CONFIG_PATH="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/2.2.busco/augustus/config/"


## run BUSCO - with arthropoda_odb9 ====================
echo "now running busco"
python /n/home12/fseixas/software/busco_v3/scripts/run_BUSCO.py \
-i $FAS \
--cpu 8 \
--out $SPP.HM.arthropoda_odb9 \
--mode genome \
--blast_single_core \
--lineage_path /n/mallet_lab/fseixas/1.projects/1.pseudo_references/2.2.busco/arthropoda_odb9/ \
--tmp_path /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/2.2.busco/tmp.$SPP/

