# create necessary directories
mkdir assemblies
mkdir hmelv25
mkdir heradem

# copy asemblies
for file in `ls /n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/0.1.HaploMerger/HaploMerger2_20180603/*.UPqscore/*ref.fa.gz`; do cp $file assemblies/; done
for file in `ls assemblies/*`; do gunzip $file; done
