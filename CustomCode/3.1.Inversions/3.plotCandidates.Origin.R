# load libraries
require(ggplot2)


## variables ==================================================
ref="hmelv25"
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
minInvSize = 100000


## 
# define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")
# define chromosomes
if (ref =="hmelv25") { bigchroms = paste0("Hmel2", sprintf("%02d", seq(1,21))) }
if (ref =="heradem") { bigchroms = paste0("Herato", sprintf("%02d", seq(1,21))) }
if (ref =="heralat") { bigchroms = paste0("Chr_", seq(1,21)) }


## read and filter ==================================================
# read candidates list
InvCandidates.BreakpointsFilter = read.table(paste0(dir,"InvCandidates.BreakpointsFilter-2-",ref,".txt"), header=T)
InvCandidates.BreakpointsFilter$species = factor(InvCandidates.BreakpointsFilter$species, levels = c(mel_clade,mid_clade,era_clade))
# filter candidates based on size
ppp1 = subset(InvCandidates.BreakpointsFilter, Bend-Bsta >= minInvSize)
ppp2 = subset(ppp1, Bend-Bsta < 3000000)
table(InvCandidates.BreakpointsFilter$species)
# ggplot(ppp2, aes(x=species, y=(Bend-Bsta)/1000)) +
#   geom_boxplot(alpha=0.25) + geom_jitter(alpha=0.5) + coord_flip() +
#   xlab("Species") + ylab("Candidate Inversion Size (Kb)")


## count how many times each inversion is supported (in all species together) ===========================
ppp2$species = ifelse(ppp2$species == "hhimfat", "hhim", as.character(ppp2$species))
ppp2$species = factor(ppp2$species, levels = c(mel_clade,mid_clade,era_clade))
#
spSummary = data.frame()
buffer = 50000
for (ch in unique(ppp2$R.scaffold)) {
  s1 = subset(ppp2, R.scaffold == ch)
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

## apply filter: #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# spSummaryFilter = subset(spSummary, count > 1 & count <= 8)
spSummaryFilter = subset(spSummary, count >= 2 & spCnt >= 2)
table(droplevels(spSummaryFilter$R.scaffold))



## arrange data so that same chromosome scaffolds come together ========================================
if (ref == "hmelv25") {
  spSummaryFilter$chrom = substr(spSummaryFilter$R.scaffold, 1, 7)
  spSummaryFilter$scaff = substr(spSummaryFilter$R.scaffold, 8, 11)
}
if (ref == "heradem") {
  spSummaryFilter$chrom = substr(spSummaryFilter$R.scaffold, 1, 8)
  spSummaryFilter$scaff = substr(spSummaryFilter$R.scaffold, 9, 10)
}
if (ref == "heralat") {
  spSummaryFilter$chrom = sapply(str_split(spSummaryFilter$R.scaffold, "_"), "[", 2)
  spSummaryFilter$scaff = sapply(str_split(spSummaryFilter$R.scaffold, "_"), "[", 3)
}

# get scaffold start and end positions relative to chromosome ////////////////////
scaffLengths = read.table(paste0("/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/",ref,".scaffoldsLength.txt"))
names(scaffLengths) = c("scaffold", "length")
if (ref == "hmelv25") { scaffLengths$chrom = substr(scaffLengths$scaffold, 1, 7) }
if (ref == "heradem") { scaffLengths$chrom = substr(scaffLengths$scaffold, 1, 8) }
if (ref == "heralat") { scaffLengths$chrom = sapply(str_split(scaffLengths$scaffold, "_"), "[", 2) }

# scaffLengths$CumulativeLength = scaffLengths$length
scaffLengths$c_sta = 0
scaffLengths$c_end = 0
# define first line cumulative coordinates
scaffLengths$c_sta[1] = 1
scaffLengths$c_end[1] = scaffLengths$length[1]
# update cumulative coordinates
for (rw in 2:nrow(scaffLengths)) {
  # in the same chromosome as previous scaffold
  if (scaffLengths$chrom[rw] == scaffLengths$chrom[rw-1]) {
    scaffLengths$c_sta[rw] = scaffLengths$c_sta[rw] + scaffLengths$c_end[rw-1] + 1
    scaffLengths$c_end[rw] = scaffLengths$c_sta[rw] + scaffLengths$length[rw]
  }
  # start of a different chromosome 
  if (scaffLengths$chrom[rw] != scaffLengths$chrom[rw-1]) {
    scaffLengths$c_sta[rw] = 1
    scaffLengths$c_end[rw] = scaffLengths$length[rw]
  }
}

## update Inversion coordinates /////
spSummaryFilter$c_Bsta = 0
spSummaryFilter$c_Bend = 0
for (rw in 1:nrow(spSummaryFilter)) {
  scaff = as.character(spSummaryFilter[rw,]$R.scaffold)
  csta = subset(scaffLengths, scaffold == scaff)$c_sta
  spSummaryFilter$c_Bsta[rw] = spSummaryFilter$Bsta[rw] + csta
  spSummaryFilter$c_Bend[rw] = spSummaryFilter$Bend[rw] + csta
}



## filter Inversions that are right next to scaffold boundaries ////////////////////
## [likely to be artifacts] ////////////////////////////////////////////////////////
spSummaryFilter2 = data.frame()
minAllowedDistance  = 50000
for (rw in 1:nrow(spSummaryFilter)) {
  # get corresponding scaffold left and right coordinates
  scaff = as.character(spSummaryFilter[rw,]$R.scaffold)
  llim = subset(scaffLengths, scaffold == scaff)$c_sta
  rlim = subset(scaffLengths, scaffold == scaff)$c_end
  # get Inversion coordinates
  ista = spSummaryFilter[rw,]$c_Bsta
  iend = spSummaryFilter[rw,]$c_Bend
  # distance to scaffold limits
  minDist = min(abs(ista-llim), abs(ista-rlim), abs(iend-llim), abs(iend-rlim))
  if (minDist >= minAllowedDistance) {spSummaryFilter2 = rbind(spSummaryFilter2, spSummaryFilter[rw,])}
}
nrow(spSummaryFilter2)
nrow(spSummaryFilter)
# write to file
write.table(file=paste0(dir,"InvCandidates.BreakpointsFilter2-2-",ref,".txt"), spSummaryFilter2, quote=F, row.names=F, sep="\t")


# ggplot(spSummaryFilter2) + 
#   geom_bar(aes(chrom)) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))


## PLOT INVERSIONS ========================================

## plot ALL chromosomes /////////////////////////////////////////////////////
# all species
ggplot(spSummaryFilter2, aes(label=Q.scaffold)) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, 
                height=0.8, fill=species)) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3)

# subset of species 
if (ref == "hmelv25") {sppFilter = subset(spSummaryFilter2, species %in% mel_clade)}
if (ref == "heradem") {sppFilter = subset(spSummaryFilter2, species %in% era_clade)}

ggplot(sppFilter, aes(label=Q.scaffold)) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, 
    width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species)) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3)

subset(sppFilter, chrom == bigchroms[15])

nrow(subset(sppFilter, Bend-Bsta > 100000))

ggplot(spSummaryFilter2, aes(label=Q.scaffold)) +  
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, 
                height=0.8, fill=species)) +
  geom_vline(data=scaffLengths, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~chrom, nrow=3)


# plot specific chromosome
chr = bigchroms[17]
spSummaryFilterSub = subset(spSummaryFilter, chrom == chr)
scaffLengthsSub = subset(scaffLengths, chrom == chr)
# plot all species
# xmin = (min(spSummaryFilterSub$Bsta)-50000)/1000000
# xmax = (max(spSummaryFilterSub$Bend)+50000)/1000000
ggplot(spSummaryFilterSub) +
  geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), col="black", alpha=0.5) +
  geom_vline(data=scaffLengthsSub, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
  geom_vline(data=scaffLengthsSub, aes(xintercept = c_end/1000000), col="grey", lty=2) +
  # xlim(xmin,xmax) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggtitle(chr) + 
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

