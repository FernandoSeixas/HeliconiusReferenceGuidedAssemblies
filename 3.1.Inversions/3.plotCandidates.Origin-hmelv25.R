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
era_clade = c("hera","hhimfat","hher","hhim","hsia","htel","hdem","hsar","hnat")
# define chromosomes
if (ref =="hmelv25") { bigchroms = paste0("Hmel2", sprintf("%02d", seq(1,21))) }
if (ref =="heradem") { bigchroms = paste0("Herato", sprintf("%02d", seq(1,21))) }
if (ref =="heralat") { bigchroms = paste0("chr", seq(1,21)) }
# limits for plots
if (ref == "hmelv25") { ylims = c(8.5,10.5) }
if (ref == "heradem") { ylims = c(8.5,10.5) }
if (ref == "heralat") { ylims = c(8.5,10.5) }

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



#################### ADD FILTERS // DIFF COUNTS ABOUT SPP WITH INV #########################
####################     ALL SCAFFOLDS SUPPORTING IT, ETC ...      #########################

ICVsFilter1 = subset(ICVs, len >= minInvSize & len <= maxInvSize)

# ICVsFilter1$inReference = NA
ICVsFilter1$nbSpecies = NA
ICVsFilter1$SameSpecies = NA
ICVsFilter1$CountAll = NA

ICVsFilter1$melClade = NA
ICVsFilter1$midClade = NA
ICVsFilter1$eraClade = NA

if (ref == "hmelv25") { refSpeciesInversions = subset(ICVsFilter1, species == "hmel") }
if (ref == "heradem") { refSpeciesInversions = subset(ICVsFilter1, species == "hera") }

# bf = 25000
# for (rw in 1:nrow(ICVsFilter1)) {
#   s1 = ICVsFilter1[rw,]
#   sp = as.character(s1$species)
#   sc = as.character(s1$R.scaffold)
#   st = s1$Bsta
#   en = s1$Bend
#   ICVsFilter1$inReference[rw] = nrow(unique(subset(refSpeciesInversions, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf)))
#   ICVsFilter1$nbSpecies[rw] = length(unique(subset(ICVsFilter1, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf)$species))
#   ICVsFilter1$SameSpecies[rw] = nrow(subset(ICVsFilter1, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species == sp))
#   ICVsFilter1$CountAll[rw] = nrow(subset(ICVsFilter1, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf))
#   
#   ICVsFilter1$melClade[rw] = length(unique(subset(ICVsFilter1, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species %in% mel_clade)$species))
#   ICVsFilter1$midClade[rw] = length(unique(subset(ICVsFilter1, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species %in% mid_clade)$species))
#   ICVsFilter1$eraClade[rw] = length(unique(subset(ICVsFilter1, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species %in% era_clade)$species))
# }

for (rw in 1:nrow(ICVsFilter1)) {
  # define variables
  s1 = ICVsFilter1[rw,]
  sp = as.character(s1$species)
  sc = as.character(s1$R.scaffold)
  st = s1$Bsta
  en = s1$Bend
  # get other overlapping inversion candidates 
  s2 = subset(ICVsFilter1, R.scaffold == sc & st < Bend & en > Bsta)
  en-st
  # s2$overlap = -1
  # s2$maxSize = -1
  # s2$percOvl = -1
  for (i in 1:nrow(s2)) {
    s2$overlap[i] = abs(min(en, s2$Bend[i])-max(st, s2$Bsta[i]))+1
    s2$maxSize[i] = max((s2$Bend[i]-s2$Bsta[i]),(en-st))+1
    s2$percOvl[i] = s2$overlap[i]/s2$maxSize[i]
  }
  s2 = subset(s2, percOvl >= 0.75)
  # count specific patterns
  ICVsFilter1$CountAll[rw] = nrow(s2)
  ICVsFilter1$nbSpecies[rw] = length(unique(s2$species))
  ICVsFilter1$SameSpecies[rw] = nrow(subset(s2, species == sp))
  ICVsFilter1$melClade[rw] = length(unique(subset(s2,species %in% mel_clade)$species))
  ICVsFilter1$midClade[rw] = length(unique(subset(s2,species %in% mid_clade)$species))
  ICVsFilter1$eraClade[rw] = length(unique(subset(s2,species %in% era_clade)$species))
}
nrow(ICVsFilter1)

#################### DEFINE SPECIES Y-POS ####################
ICVsFilter1$species_ypos = 0
covTable$species_ypos = 0
# define ypos
for (i in 1:length(levels(ICVsFilter1$species))) {
  spp = levels(ICVsFilter1$species)[i]
  ICVsFilter1$species_ypos = ifelse(ICVsFilter1$species == spp, i, ICVsFilter1$species_ypos)
  covTable$species_ypos = ifelse(covTable$species == spp, i, covTable$species_ypos)
}
# ylims
if (ref == "hmelv25") { ylims = c(8.5, 10.5) }
if (ref == "heradem") { ylims = c(8.5, 10.5) }

# mmm = subset(ICVsFilter1, chrom == "Herato21")
# mmm[order(mmm$Bsta),]



######################### PLOT #########################
# ICVsFilter2 = subset(ICVsFilter1, nbSpecies >= 2)

# # filter by groups
ICVsFilter2 = subset(ICVsFilter1, melClade <= 3 & eraClade >= 3)
# ICVsFilter2 = subset(ICVsFilter1, melClade >= 3 & eraClade <= 2)


ICVsFilter2 = subset(ICVsFilter1, melClade >= 2 & eraClade == 0 & midClade == 0 & len < maxInvSize)
# ICVsFilter2 = subset(ICVsFilter1, melClade >= 2 & eraClade == 0 & len > 50000)
# ICVsFilter2 = subset(ICVsFilter1, melClade >= 2  & nbSpecies < 7 | eraClade >= 2 & nbSpecies < 7)

# ICVsFilter2 = subset(ICVsFilter1, SameSpecies >= 2 & melClade >= 2 & eraClade <= 0)
# ICVsFilter2 = subset(ICVsFilter1, SameSpecies >= 2 & eraClade >= 2)

# ICVsFilter2 = subset(ICVsFilter1, nbSpecies >= 2 & nbSpecies < 7)
# ICVsFilter3 = subset(ICVsFilter2, melClade >= 2 & eraClade == 0)

# plot to check distribution of candidates =============
ggplot(ICVsFilter1) +
# ggplot(ICVsFilter2) +
  # candidate inversions
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.8) +
  # add coverage
  # geom_point(data=covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, colour=covOutlier, group=species)) +
  # scale_color_manual(values = c(alpha("black",0.0),alpha("red",0.1))) +
  # add reference scaffold limits
  geom_vline(data=scaffLengths, aes(xintercept=c_sta/1000000), lty=2, col="grey") +
  geom_vline(data=scaffLengths, aes(xintercept=c_end/1000000), lty=2, col="grey") +
  # species limits
  geom_hline(yintercept = ylims, alpha=0.5) +
  # add hmelv25/heradem/heralat
  # geom_gene_arrow(data=ref2ref_1_Merge, aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=18), fill=ifelse(ref2ref_1_Merge$strand == "-", "black", "#f5f5f5"), lwd=0.0) +
  # geom_gene_arrow(data=ref2ref_2_Merge, aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=19), fill=ifelse(ref2ref_2_Merge$strand == "-", "black", "#f5f5f5"), lwd=0.0) +
  # visual
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3, scales="free_x")



## plot specific chromosome ============================
# subset data
chr = bigchroms[6]
# ICVsFilter3 = subset(ICVsFilter1, chrom == chr)
# ICVsFilter3 = subset(ICVsFilter1, chrom == chr & len < maxInvSize)
# ICVsFilter3 = subset(ICVsFilter1, chrom == chr & nbSpecies >= 2 & len < 3000000) 
sub_covTable=subset(covTable, chrom == chr)
sub_scaffLen=subset(scaffLengths, chrom == chr)
# sub_ref2ref_1 = subset(ref2ref_1_Merge, chrom == chr & strand == "-" & abs(R.end-R.start) > 5000)
# sub_ref2ref_2 = subset(ref2ref_2_Merge, chrom == chr & strand == "-" & abs(R.end-R.start) > 5000)
sub_ref2ref_1 = subset(ref2ref_1_Merge, chrom == chr)
sub_ref2ref_2 = subset(ref2ref_2_Merge, chrom == chr)
# plot
ggplot(ICVsFilter3) +
  # candidate inversions
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  # add hmelv25/heradem/heralat
  geom_gene_arrow(data=sub_ref2ref_1, aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=18, color=Q.scaffold), fill=ifelse(sub_ref2ref_1$strand == "-", "black", "#f5f5f5"), lwd=0.6) +
  geom_gene_arrow(data=sub_ref2ref_2, aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=19, color=Q.scaffold), fill=ifelse(sub_ref2ref_2$strand == "-", "black", "#f5f5f5"), lwd=0.6) +
  # add ref scaffold limits
  geom_vline(xintercept = sub_scaffLen$c_sta/1000000, col="grey", lty=2) +
  geom_vline(xintercept = sub_scaffLen$c_end/1000000, col="grey", lty=2) +
  # add coverage
  # geom_point(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, colour=covOutlier, group=species)) +
  # scale_color_manual(values = c(alpha("black",0.0),alpha("red",0.1))) +
  # geom_line(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos+relcov/max(sub_covTable$relcov), group=species), alpha=0.5) +
  geom_line(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos+relcov/10, group=species), alpha=0.5) +
  # geom_line(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos+relcov, group=species), alpha=0.3) +
  # species limits
  geom_hline(yintercept = ylims) +
  # geom_vline(xintercept = c(2553006,9375339)/1000000) +
  # xlim(2.5,3.5) +
  # ylim(0,20) +
  ylab("") +
  xlab("Chromosome Position (Mb)") + ggtitle(chr)







######################### summarize all information #########################

##  ==================================================
GlobalSummary = data.frame()
for (rscaff in unique(ICVsFilter1$R.scaffold)) {
  s1 = subset(ICVsFilter1, R.scaffold == rscaff)
  s1.order = s1[order(s1$Bsta),]
  s1.order$count = 1
  s1.merge = s1.order[1,]
  s1.merge$melClade = 0
  s1.merge$midClade = 0
  s1.merge$eraClade = 0
  #
  # if (nrow)
  if (nrow(s1.order) > 1) {
    for (rw in 2:nrow(s1.order)) {
      Mrow = nrow(s1.merge)
      x1 = s1.merge$Bsta[Mrow] ; x2 = s1.merge$Bend[Mrow]
      y1 = s1.order$Bsta[rw] ;   y2 = s1.order$Bend[rw]
      overlap = min(x2, y2)-max(x1, y1)+1 # length of overlap
      maxSize = max((x2-x1),(y2-y1))
      percOvl = overlap/maxSize
      # same inversion
      if (percOvl >= 0.75) {
        # update coordinates and count
        s1.merge$Bend[Mrow] = max(x1,y2)
        s1.merge$c_Bend[Mrow] = max(s1.order$c_Bend[rw],s1.merge$c_Bend[Mrow])
        s1.merge$count[Mrow] = s1.merge$count[Mrow] + 1
        if (s1.order$species[rw] %in% mel_clade) { s1.merge$melClade[Mrow] = s1.merge$melClade[Mrow] + 1 }
        if (s1.order$species[rw] %in% mid_clade) { s1.merge$midClade[Mrow] = s1.merge$midClade[Mrow] + 1 }
        if (s1.order$species[rw] %in% era_clade) { s1.merge$eraClade[Mrow] = s1.merge$eraClade[Mrow] + 1 }
      }
      # diff inversion
      if (percOvl < 0.75) { s1.merge = rbind(s1.merge, s1.order[rw,])}
    } 
  }
  GlobalSummary = rbind(GlobalSummary, s1.merge)
}

plot(table(GlobalSummary$count))

nrow(subset(GlobalSummary, count == 1))
nrow(subset(GlobalSummary, count >= 2))
nrow(subset(GlobalSummary, count <= 20))

subset(GlobalSummary, count > 20)

mmm = subset(GlobalSummary, melClade >= 3 & eraClade == 0 & midClade == 0)

## plot ======================================================================
ggplot(ICVsFilter1) +
  # candidate inversions
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, width=(c_Bend-c_Bsta+200000)/1000000, height=0.8, fill=species)) +
  # add inversions summary
  geom_tile(data=GlobalSummary, aes(x=(c_Bsta+c_Bend)/2000000, y=0, width=(c_Bend-c_Bsta)/1000000, height=0.8), alpha=0.5) +
  # geom_tile(data=subset(GlobalSummary, count >= 3), aes(x=(c_Bsta+c_Bend)/2000000, y=0, width=(c_Bend-c_Bsta)/1000000, height=0.8), alpha=0.5) +
  # add reference scaffold limits
  geom_vline(data=scaffLengths, aes(xintercept=c_sta/1000000), lty=2, col="grey") +
  geom_vline(data=scaffLengths, aes(xintercept=c_end/1000000), lty=2, col="grey") +
  # species limits
  geom_hline(yintercept = ylims, alpha=0.5) +
  # visual
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3, scales="free_x")



#####
read.table(dir)
