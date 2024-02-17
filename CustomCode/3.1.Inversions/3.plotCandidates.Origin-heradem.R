# load libraries
require(ggplot2)
require(gggenes)

## variables ==================================================
ref="heradem"
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
# define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")
# define chromosomes
if (ref =="hmelv25") { bigchroms = paste0("Hmel2", sprintf("%02d", seq(1,21))) }
if (ref =="heradem") { bigchroms = paste0("Herato", sprintf("%02d", seq(1,21))) }
if (ref =="heralat") { bigchroms = paste0("Chr_", seq(1,21)) }
# limits for plots
if (ref == "hmelv25") { ylims = c(8.5,10.5) }
if (ref == "heradem") { ylims = c(8.5,10.5) }
if (ref == "heralat") { ylims = c(8.5,10.5) }

# read data
ICVs = read.table(paste0(dir,"InvCandidates.Breakpoints-2-",ref,".txt"), header=T)
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
hmelv25_2_heradem_Minimap2Merge = read.table(paste0(dir,"hmelv25_2_heradem_Minimap2Merge.txt"), header=T)
heralat_2_heradem_Minimap2Merge = read.table(paste0(dir,"heralat_2_heradem_Minimap2Merge.txt"), header=T)



#################### ADD FILTERS // DIFF COUNTS ABOUT SPP WITH INV #########################
####################     ALL SCAFFOLDS SUPPORTING IT, ETC ...      #########################
ICVs$inReference = NA
ICVs$nbSpecies = NA
ICVs$SameSpecies = NA
ICVs$CountAll = NA

ICVs$melClade = NA
ICVs$midClade = NA
ICVs$eraClade = NA

if (ref == "hmelv25") { refSpeciesInversions = subset(ICVs, species == "hmel") }
if (ref == "heradem") { refSpeciesInversions = subset(ICVs, species == "hera") }

bf = 25000
for (rw in 1:nrow(ICVs)) {
  s1 = ICVs[rw,]
  sp = as.character(s1$species)
  sc = as.character(s1$R.scaffold)
  st = s1$Bsta
  en = s1$Bend
  ICVs$inReference[rw] = nrow(unique(subset(refSpeciesInversions, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf)))
  ICVs$nbSpecies[rw] = length(unique(subset(ICVs, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf)$species))
  ICVs$SameSpecies[rw] = nrow(subset(ICVs, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species == sp))
  ICVs$CountAll[rw] = nrow(subset(ICVs, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf))
  
  ICVs$melClade[rw] = length(unique(subset(ICVs, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species %in% mel_clade)$species))
  ICVs$midClade[rw] = length(unique(subset(ICVs, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species %in% mid_clade)$species))
  ICVs$eraClade[rw] = length(unique(subset(ICVs, R.scaffold == sc & abs(Bsta-st) < bf & abs(Bend-en) < bf & species %in% era_clade)$species))
}



#################### DEFINE SPECIES Y-POS ####################
ICVs$species_ypos = 0
covTable$species_ypos = 0
# define ypos
for (i in 1:length(levels(ICVs$species))) {
  spp = levels(ICVs$species)[i]
  ICVs$species_ypos = ifelse(ICVs$species == spp, i, ICVs$species_ypos)
  covTable$species_ypos = ifelse(covTable$species == spp, i, covTable$species_ypos)
}



######################### PLOT #########################
ICVsFilter1 = subset(ICVs, len >= minInvSize & len < maxInvSize)

table(ICVsFilter$melClade)

length(mel_clade)
length(era_clade)

# ICVsFilter = subset(ICVsFilter1, melClade >= 4 & eraClade < 2)
ICVsFilter = subset(ICVsFilter1, melClade <= 2 & eraClade >= 2)

# plot to check distribution of candidates =============
ggplot(ICVsFilter) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  # add ref scaffold limits
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=scaffLengths, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  # add coverage
  # geom_point(data=covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, colour=covOutlier, group=species)) +
  # scale_color_manual(values = c(alpha("black",0.0),alpha("red",0.1))) +
  # visual
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3, scales="free_x")

ggplot(ICVs) + geom_bar(aes(species))

## plot specific chromosome ============================
# subset data
sc = "Herato13"
sub_ICVs = subset(ICVsFilter, chrom == sc & len >= minInvSize & len < maxInvSize )
sub_covTable = subset(covTable, chrom == sc)
sub_scaffLen = subset(scaffLengths, chrom == sc)
# sub_Repeats = subset(Repeats, chrom == sc)
# plot
ggplot() +  
  geom_tile(data=sub_ICVs, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  # add ref scaffold limits
  geom_vline(data=sub_scaffLen, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=sub_scaffLen, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  # add coverage
  geom_point(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, colour=covOutlier, group=species)) +
  geom_line(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos+(relcov/max(sub_covTable$relcov)), group=species)) +  
  scale_color_manual(values = c(alpha("black",0.0),alpha("red",0.5))) +
  # geom_tile(data=sub_Repeats, aes(x=(c_Bsta+c_Bend)/2000000, y=-1, width=(c_Bend-c_Bsta)/1000000, height=0.8)) +
  ggtitle(sc)

sub_ICVs[order(sub_ICVs$c_Bsta),]



######################### APPLY FILTERS #########################

minInvSize = 100000
maxInvSize = 3000000
ICVs = subset(ICVs, len >= minInvSize & len < maxInvSize)
head(ICVs)

## filter
ICVsFilter = subset(ICVs, inReference == 0)
ICVsFilter = subset(ICVs, SameSpecies >= 2 | CountAll >= 2)
table(ICVs$inReference)
table(ICVs$nbSpecies)
table(ICVs$SameSpecies)
table(ICVs$CountAll)

# plot all chromosomes
ICVsFilter = subset(ICVs, stCov >= 3 | enCov >= 3)

ggplot(ICVsFilter) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=scaffLengths, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  # geom_point(data=covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, colour=covOutlier, group=species)) +
  # scale_color_manual(values = c(alpha("black",0.0),alpha("red",0.1))) +
  # geom_line(data=covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos+relcov/max(covTable$relcov), group=species), alpha=0.5) +
  facet_wrap(~chrom, nrow=3, scales="free_x")


# filter 
chr = "Herato17"
ICVsFilter = subset(ICVs, chrom == chr)
sub_covTable=subset(covTable, chrom == chr)
sub_scaffLen=subset(scaffLengths, chrom == chr)
sub_hmelv25 = subset(hmelv25_2_heradem_Minimap2Merge, chrom == chr & strand == "-" & abs(R.end-R.start) > 5000)
sub_heralat = subset(heralat_2_heradem_Minimap2Merge, chrom == chr & strand == "-" & abs(R.end-R.start) > 5000)
# plot
ggplot(ICVsFilter) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  # add chromosome breakpoints
  geom_vline(xintercept = sub_scaffLen$c_sta/1000000, col="grey", lty=2) +
  geom_vline(xintercept = sub_scaffLen$c_end/1000000, col="grey", lty=2) +
  # add coverage
  geom_point(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos, colour=covOutlier, group=species)) +
  scale_color_manual(values = c(alpha("black",0.0),alpha("red",0.1))) +
  geom_line(data=sub_covTable, aes(x=(c_Bsta+c_Bend)/2000000, y=species_ypos+relcov/max(sub_covTable$relcov), group=species), alpha=0.5) +
  # add hmelv25/heralat
  geom_gene_arrow(data=sub_hmelv25, aes(xmin = R.start/1000000, xmax = R.end/1000000, y=1)) +
  geom_gene_arrow(data=sub_heralat, aes(xmin = R.start/1000000, xmax = R.end/1000000, y=18)) +
  #
  facet_wrap(~chrom, nrow=3, scales="free_x")

#
ICVsFilter = subset(ICVs, SameSpecies == 4)

ggplot(ICVsFilter) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=Q.scaffold, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6)




refSpeciesCovTable = subset(covTable, species == "hmel")

ggplot(refSpeciesInversions) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=1, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  geom_line(data=refSpeciesCovTable, aes(x=(c_Bsta+c_Bend)/2000000, y=relcov)) +
  facet_wrap(~chrom, nrow=3, scales="free_x")
tail(refSpeciesCovTable)



################### count how many times each inversion is supported ###################
################### (in all species together) ##########################################
#
spSummary = data.frame()
buffer = 50000
for (ch in unique(ICVs$R.scaffold)) {
  s1 = subset(ICVs, R.scaffold == ch)
  s1 = s1[order(s1$Bsta),]
  s1$count = 0
  s1$spCnt = 0
  if (nrow(s1) >= 1) { 
    for (rw in 1:nrow(s1)) {
      st = s1$Bsta[rw]
      en = s1$Bend[rw]
      ss = subset(s1, Bsta > st-buffer & Bsta < st+buffer)
      ss = subset(ss, Bend > en-buffer & Bend < en+buffer)
      s1$count[rw] = nrow(ss)
      s1$spCnt[rw] = length(unique(ss$species))
    }
  }
  spSummary = rbind(spSummary, s1)
}
# apply filter: chose only inversions found in more than one species 
# spSummaryFilter = subset(spSummary, count > 1 & count <= 8)
spSummaryFilter = subset(spSummary, count >= 1 | spCnt >= 1)
spSummaryFilter = subset(spSummary, count >= 2 | spCnt >= 2)
table(droplevels(spSummaryFilter$R.scaffold))

nrow(spSummaryFilter)

## count scaffolds that support inversion per species ======================================
maxBreakDist = 25000
for (rw in 1:nrow(ICVs)) {
  s1 = ICVs[rw,]
  sppInv = s1$species
  chrInv = s1$R.scaffold
  staInv = s1$Bsta
  endInv = s1$Bend
  ICVs$SameSpecies[rw] = nrow(subset(ICVs, species == sppInv & R.scaffold == chrInv & abs(staInv-ICVs$Bsta) < maxBreakDist & abs(endInv-ICVs$Bend) < maxBreakDist))
}
# filter
BreakSupportFilter = subset(ICVs, SameSpecies == 2 | SameSpecies == 4)
nrow(BreakSupportFilter)


## PLOT INVERSIONS ==========================================================
#

if (ref == "hmelv25") { spSummaryFilter$chrom =  factor(spSummaryFilter$chrom, levels = paste0("Hmel2", sprintf("%02d", seq(1,21)))) }
if (ref == "heradem") { spSummaryFilter$chrom =  factor(spSummaryFilter$chrom, levels = paste0("Herato", sprintf("%02d", seq(1,21)))) }
if (ref == "heralat") { spSummaryFilter$chrom =  factor(spSummaryFilter$chrom, levels = paste0("chr", seq(1,21))) }

if (ref == "hmelv25") { BreakSupportFilter$chrom =  factor(BreakSupportFilter$chrom, levels = paste0("Hmel2", sprintf("%02d", seq(1,21)))) }
if (ref == "heradem") { BreakSupportFilter$chrom =  factor(BreakSupportFilter$chrom, levels = paste0("Herato", sprintf("%02d", seq(1,21)))) }
if (ref == "heralat") { BreakSupportFilter$chrom =  factor(BreakSupportFilter$chrom, levels = paste0("chr", seq(1,21))) }

## all species ==================================================
# all chromosomes
ggplot(spSummaryFilter, aes(label=Q.scaffold)) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), alpha=0.6) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=scaffLengths, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  geom_hline(yintercept = ylims, col="grey") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3, scales="free_x")

ggplot(BreakSupportFilter, aes(label=Q.scaffold)) +
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species)) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=scaffLengths, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3, scales="free_x")

ggplot(BreakSupportFilter, aes(label=Q.scaffold)) +
  # geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species)) +
  # geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  # geom_vline(data=scaffLengths, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3, scales="free_x")


ggplot(BreakSupportFilter) + geom_bar(aes(species)) + scale_y_continuous(breaks=seq(0,20,1))

a = as.data.frame(table(BreakSupportFilter$species))
b = aggregate(minimaps$Q.length, list(minimaps$species), median)
names(a) = c("species", "inversions")
names(b) = c("species", "medianScaffSize")
c = merge(a, b, by = "species")
ggplot(c) + geom_point(aes(x=inversions, y=medianScaffSize/1000))




# plot specific chromosome
chr = bigchroms[17]
spSummaryFilterSub = subset(spSummaryFilter, chrom == chr)
scaffLengthsSub = subset(scaffLengths, chrom == chr)
ggplot(spSummaryFilterSub) +
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), col="black", alpha=0.5) +
  geom_vline(data=scaffLengthsSub, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=scaffLengthsSub, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  scale_x_continuous(breaks = seq(0,25,1)) +
  # xlim(xmin,xmax) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggtitle(chr) + 
  xlab("Chromosome Position(Mb)") + ylab("Species")


subset(spSummaryFilterSub, species == "htel")
subset(spSummaryFilterSub, species == "hsia")
subset(InvCandidates.Pass1, species == "hpar" & Q.scaffold == "Sc0000663")
subset(InvCandidates.Pass1, species == "hele" & Q.scaffold == "Sc0001747")

mmm = subset(minimaps, Q.scaffold == "xfSc0002464" & species == "hsar")
mmm = subset(minimaps, Q.scaffold == "xfSc0001066" & species == "hdem")
mmm = subset(minimaps, Q.scaffold == "Sc0000432" & species == "hdem")
mmm = subset(minimaps, Q.scaffold == "xfSc0004874" & species == "htel")
mmm = subset(minimaps, Q.scaffold == "xfSc0000508" & species == "htel")
mmm = subset(minimaps, Q.scaffold == "xfSc0001133" & species == "hsia")
mmm = subset(minimaps, Q.scaffold == "xfSc0001618" & species == "hsia")
mmm[order(mmm$Q.start),]


# subset of species 
if (ref == "hmelv25") {sppFilter = subset(spSummaryFilter, species %in% mel_clade)}
if (ref == "heradem") {sppFilter = subset(spSummaryFilter, species %in% era_clade)}

ggplot(sppFilter, aes(label=Q.scaffold)) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, 
                width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species)) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3)

spSummaryFilter2 = subset(sppFilter, chrom == bigchroms[21])
ggplot(spSummaryFilter2, aes(label=Q.scaffold)) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, 
                height=0.8, fill=species)) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species")



spSummaryFilterSub[order(spSummaryFilterSub$Bsta),]


# # sppFilter = subset(spSummaryFilterSub, species %in% c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar"))
# xmin = (min(sppFilter$Bsta)-50000)/1000000
# xmax = (max(sppFilter$Bend)+50000)/1000000
# ggplot(sppFilter, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
#   geom_tile(col="black", alpha=0.5) +
#   geom_vline(data=refBreaksSub, xintercept = refBreaksSub$cumSta/1000000) +
#   xlim(xmin, xmax) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
#   ggtitle(chr) + xlab("Chromosome Position(Mb)") + ylab("Species")
# 
# sppFilter[order(sppFilter$Bsta),]



## summarize all information ==================================================
spSummaryFilter = subset(spSummary, spCnt >= 2)

buffer=25000
AllSummary = data.frame()
for (scaff in unique(spSummaryFilter$R.scaffold)) {
  s1 = subset(spSummaryFilter, R.scaffold == scaff)
  s1 = s1[order(s1$Bsta),]
  s1$species = as.character(s1$species)
  s1$count = 1
  s1$spCnt = 1
  s2 = s1[1,]
  for (rw in 2:nrow(s1)) {
    if (s1$Bsta[rw]-s2$Bsta[nrow(s2)] < buffer & s1$Bend[rw]-s2$Bend[nrow(s2)] < buffer) {
      s2$Bsta[nrow(s2)] = min(s1$Bsta[rw],s2$Bsta[nrow(s2)])
      s2$Bend[nrow(s2)] = max(s1$Bend[rw],s2$Bend[nrow(s2)])
      s2$count[nrow(s2)] = s2$count[nrow(s2)] + 1
      s2$spCnt[nrow(s2)] = s2$spCnt[nrow(s2)] + 1
      if (s2$species[nrow(s2)] != s1$species[rw]) {
        s2$species[nrow(s2)] = paste0(s2$species[nrow(s2)],",",s1$species[rw])
      }
    }
    if (s1$Bsta[rw]-s2$Bsta[nrow(s2)] >= buffer | s1$Bend[rw]-s2$Bend[nrow(s2)] >= buffer) {
      s2 = rbind(s2, s1[rw,])
    }
  }
  for (i in 1:nrow(s2)) {
    s2$species[i] = paste(unique(strsplit(s2$species[i], ",")[[1]]), collapse = ",")
    s2$spCnt[i] = length(unique(strsplit(s2$species[i], ",")[[1]]))
  }
  AllSummary = rbind(AllSummary, s2)
}

AllSummaryFilter = subset(AllSummary, spCnt >= 2)
subset(AllSummary, R.scaffold == "Herato1701" & spCnt >= 2)


ggplot(AllSummaryFilter, aes(x = (Bsta+Bend)/2000000, y=1, width=(Bend-Bsta)/1000000, height=0.8)) +
  geom_tile(col="black", alpha=0.5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~R.scaffold) + xlab("Chromosome Position(Mb)") + ylab("Species")

