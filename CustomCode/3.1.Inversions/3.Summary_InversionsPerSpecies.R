# load libraries
require(ggplot2)
library(gggenes)
require(stringr)

## variables ==================================================
minInvSize = 50000
maxInvSize = 2000000

ref="heradem"
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/"
# define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")
# define chromosomes
if (ref =="hmelv25") { bigchroms = paste0("Hmel2", sprintf("%02d", seq(1,21))) }
if (ref =="heradem") { bigchroms = paste0("Herato", sprintf("%02d", seq(1,21))) }
if (ref =="heralat") { bigchroms = paste0("chr", seq(1,21)) }

# read data
ICVs = read.table(paste0(dir,"InvCandidates.Breakpoints-2-",ref,".txt"), header=T)
nrow(ICVs)

# update species
# ICVs$species = ifelse(ICVs$species == "hhimfat", "hhim", as.character(ICVs$species))
ICVs$species = factor(ICVs$species, levels = c(mel_clade,mid_clade,era_clade))
ICVs$len = ICVs$Bend - ICVs$Bsta + 1

# ref chromosome coordinates
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)

# coverage data 
covTable = read.table(paste0(dir,"coverageTable.",ref,".txt"), header=T)
covTable$covOutlier = ifelse(covTable$relcov > 2, "red", "black")

# read ref-2-ref mappings
if (ref == "hmelv25") {
  ref2ref_1_Merge = read.table(paste0(dir,"heradem_2_hmelv25_Minimap2Merge.txt"), header=T)
  ref2ref_2_Merge = read.table(paste0(dir,"heralat_2_hmelv25_Minimap2Merge.txt"), header=T)
}
if (ref == "heradem") {
  ref2ref_1_Merge = read.table(paste0(dir,"hmelv25_2_heradem_Minimap2Merge.txt"), header=T)
  ref2ref_2_Merge = read.table(paste0(dir,"heralat_2_heradem_Minimap2Merge.txt"), header=T)
}


#
ICVsFilter1 = subset(ICVs, len >= minInvSize & len <= maxInvSize)
nrow(ICVsFilter1)

######################### INVERSIONS PER SPECIES #########################
mrgInversions = data.frame()
# go through each species
for (spp in unique(ICVsFilter1$species)) {
  s1 = subset(ICVsFilter1, species == spp)
  s1$Q.scaffold = as.character(s1$Q.scaffold)
  s1$R.scaffold = as.character(s1$R.scaffold)
  # go chrom by chromosome in the reference
  for (chrom in sort(unique(s1$R.scaffold)) ) {
    s2 = s1[s1$R.scaffold == chrom,]
    s2 = s2[order(s2$Bsta),]
    #
    mrg = s2[1,]
    if (nrow(s2) > 1) {
      # go through each w2rap scaffold
      for (rw in 2:nrow(s2)) {
        st1 = mrg$Bsta[nrow(mrg)]; en1 = mrg$Bend[nrow(mrg)]
        st2 = s2$Bsta[rw]; en2 = s2$Bend[rw]
        # overlap
        maxSize = max((en1-st1),(en2-st2))+1
        overlap = min(en1, en2)-max(st1, st2)+1
        percOvl = overlap/maxSize
        # same inversion
        if (percOvl >= 0.75) {
          stmin = min(st1, st2)
          enmax = max(en1, en2)
          mrg$Bsta[nrow(mrg)] = stmin
          mrg$Bend[nrow(mrg)] = enmax
        }
        # not same inversion
        if (percOvl < 0.75) {
          mrg = rbind(mrg, s2[rw,])
        }
      }
    }
    # add to mrgInverions table
    mrgInversions = rbind(mrgInversions, mrg) 
  }
}

nrow(mrgInversions)

# plot
ggplot(mrgInversions) + geom_bar(aes(chrom)) + coord_flip()
ggplot(mrgInversions) + geom_bar(aes(species)) + coord_flip()


table(mrgInversions$species)
min(table(mrgInversions$species))
max(table(mrgInversions$species))
