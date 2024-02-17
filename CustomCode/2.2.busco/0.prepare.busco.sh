## load modules
module load augustus/3.3-fasrc02


## prepare augustus configuration
mkdir /n/scratchlfs02/mallet_lab/fseixas/1.pseudo_references/2.2.busco/augustus
cp -r /n/helmod/apps/centos7/Core/augustus/3.3-fasrc02/config/ /n/scratchlfs02/mallet_lab/fseixas/1.pseudo_references/2.2.busco/augustus
export AUGUSTUS_CONFIG_PATH="/n/scratchlfs02/mallet_lab/fseixas/1.pseudo_references/2.2.busco/augustus/config"


## create appropriate links to data 
mkdir 0.assemblies
ln -s /n/home12/fseixas/software/busco_v3/datasets/arthropoda_odb9/ arthropoda_odb9
