## script to launch busco

# variables
ref=$1

# create dir in which to run and launch busco
mkdir $ref;
cd $ref


## Run BUSCO
for spp in `cat ~/code/heliconius_seixas/1.pseudo_references/spp.list`; do
    sbatch ~/code/heliconius_seixas/1.pseudo_references/2.2.busco/busco.HM.slurm $spp
done


## Summary of analyses
for spp in `cat ~/code/heliconius_seixas/1.pseudo_references/spp.list`; do
    printf "%s\t" $spp;
    cat run_$spp.HM.arthropoda_odb9/short_summary_$spp.HM.arthropoda_odb9.txt | grep "(S)" | awk '{printf "%s\t", $1}';  
    cat run_$spp.HM.arthropoda_odb9/short_summary_$spp.HM.arthropoda_odb9.txt | grep "(D)" | awk '{printf "%s\t", $1}';  
    cat run_$spp.HM.arthropoda_odb9/short_summary_$spp.HM.arthropoda_odb9.txt | grep "(F)" | awk '{printf "%s\t", $1}';  
    cat run_$spp.HM.arthropoda_odb9/short_summary_$spp.HM.arthropoda_odb9.txt | grep "(M)" | awk '{printf "%s\n", $1}';
done