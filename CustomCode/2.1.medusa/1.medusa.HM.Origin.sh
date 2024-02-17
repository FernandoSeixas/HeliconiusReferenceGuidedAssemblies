# variables 
spp=$1
ref=$2

# create directorie
mkdir $spp.HaploMerger
cd $spp.HaploMerger

# copy medusa default script to current folder
cp /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.1.medusa/medusa.HM.Origin.slurm .

# launch medusa
sbatch medusa.HM.Origin.slurm $spp $ref
cd ..
