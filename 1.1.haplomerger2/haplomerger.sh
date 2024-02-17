## variables
SPP=$1

## load modules and add dirs to PATH
module load perl/5.26.1-fasrc01
module load lastz/1.04.00-fasrc01
module load ucsc/20150820-fasrc01
PATH=~/software/assemblies/HaploMerger2_20180603/chainNet_jksrc20100603_centOS6/:~/software/assemblies/HaploMerger2_20180603/gapCloser_v1.12/:~/software/assemblies/HaploMerger2_20180603/SSPACE-STANDARD-3.0_linux-x86_64/:$PATH
echo $PATH


## prepare input file ===================================================
# create directory and copy necessary files
rsync -av --progress ~/software/assemblies/HaploMerger2_20180603/testProject/ ${SPP}\_DISCOVAR.UPqscore
cd ${SPP}\_DISCOVAR.UPqscore
# copy assembly to analyze
cp /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/DISCOCAR_w2rap/${SPP}.DISCOVAR.scaffolds.fa .
# softmask assembly
~/software/assemblies/HaploMerger2_20180603/winMasker/windowmasker -checkdup true -mk_counts -in ${SPP}.DISCOVAR.scaffolds.fa -out ${SPP}.DISCOVAR.masking_library.ustat -mem 6500
~/software/assemblies/HaploMerger2_20180603/winMasker/windowmasker -ustat ${SPP}.DISCOVAR.masking_library.ustat -in ${SPP}.DISCOVAR.scaffolds.fa -out ${SPP}.DISCOVAR.softmask.fasta -outfmt fasta -dust true
# DnaPolishing
cat ${SPP}.DISCOVAR.softmask.fasta | perl ~/software/assemblies/HaploMerger2_20180603/bin/faDnaPolishing.pl --legalizing --maskShortPortion=1 --noLeadingN --removeShortSeq=1 | gzip > ${SPP}\_DISCOVAR_cleaned.fa.gz


## determine Score Matrix for lastz ====================================
# load modules
module load lastz/1.04.00-fasrc01
module load perl/5.10.1-fasrc05
module load perl-modules/5.26.1-fasrc01
## convert assembly to single line fasta
zcat ${SPP}\_DISCOVAR_cleaned.fa.gz | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > ${SPP}\_DISCOVAR_cleaned.sline.fa
# see adequate size to split assembly
grep -v '>' ${SPP}\_DISCOVAR_cleaned.sline.fa | awk '{if (length($1) > 150000) sum=sum+length($1); tot=tot+length($1)} END {print sum/tot*100}'
# split assembly into big and small contigs
bash ~/code/heliconius_seixas/1.pseudo_references/2.1.medusa/filter_contigs.sh -a ${SPP}\_DISCOVAR_cleaned.sline.fa -m 150000 -b ${SPP}.part1.fa -s ${SPP}.part2.fa
# determine lastz score
#perl /n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/0.1.HaploMerger/HaploMerger2_20180603/bin/lastz_D_Wrapper.pl --target=${SPP}.part1.fa --query=${SPP}.part2.fa --identity=90
perl ~/software/assemblies/HaploMerger2_20180603/bin/lastz_D_Wrapper.pl --target=${SPP}.part1.fa --query=${SPP}.part2.fa --identity=90
# modify score matrix
tail -n 5 ${SPP}.*.q > scoreMatrix.q


## run HaploMerger2 ====================================================
# change input assembly name in files
sed -i "s/AssemblyName/${SPP}\_DISCOVAR_cleaned/g" "run_all.batch"
# run
sbatch runHaplomerger.slurm ${SPP}
# go to main folder 
cd ..
