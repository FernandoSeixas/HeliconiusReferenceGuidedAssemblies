#!/bin/bash
#SBATCH -n 8   # Number of cores requested
#SBATCH -N 1    # Ensure that all cores are on one machine
#SBATCH -t 7200 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=24000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o medusa_%j.out # Standard out goes to this file
#SBATCH -e medusa_%j.err # Standard err goes to this filehostname


## variables ============================================================
SPP=$1
REFNAM=$2
ASSEMBLY="../../assemblies/"$SPP"_DISCOVAR_cleaned_A_ref.sline.fa"
OUTPUT="$SPP-2-$REFNAM.HM.100gap.fasta"
REFFAS="/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/Origin/$REFNAM/"

## load modules 
module load Anaconda/5.0.1-fasrc02
module load mummer/4.0.0beta2-fasrc01
module load jdk/1.8.0_172-fasrc01
module load python/2.7.14-fasrc01

## running MEDUSA =======================================================
# activate medusa environment
source activate medusa_env

# copy assembly
cp $ASSEMBLY . 
INPUT=${ASSEMBLY/..\/..\/assemblies\//}

# run medusa
java -jar /n/home12/fseixas/software/medusa/medusa.jar -v \
-scriptPath ~/software/medusa/medusa_scripts \
-o $OUTPUT \
-f $REFFAS \
-i $INPUT

# abandon custom environment
source deactivate

