# load libraries
require(ggplot2)

## variables ==================================================
ref="heralat"
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
# define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("heralat", "heradem", "hera","hhimfat","hhim","hsia","htel","hdem","hsar")


## filter based on common candidates ========================================================================
InvCandidates.CoverageFilter = read.table(paste0(dir,"InvCandidates.CoverageFilter-2-",ref,".txt"), header=T)
# filter candidates based on size
ppp1 = subset(InvCandidates.CoverageFilter, Bend-Bsta >= 50000)
ppp2 = subset(ppp1, Bend-Bsta < 3000000)

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
# spSummaryFilter = subset(spSummary, count > 1 & count <= 8)
spSummaryFilter = subset(spSummary, count >= 2 & spCnt >= 2)
table(droplevels(spSummaryFilter$R.scaffold))

# plot all chromosomes
ggplot(spSummaryFilter, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~R.scaffold, nrow=4)

# plot one chromosome
if (ref == "hmelv25") {
  bigchroms = as.character(read.table("/n/mallet_lab/fseixas/1.projects/2.elevatus_pardalinus/0.data/hmelv25.scaffolds_big.txt")[1:21,1])
  spSummaryFilter$R.scaffold = factor(spSummaryFilter$R.scaffold, levels=bigchroms)
  spSummaryFilterSub = subset(spSummaryFilter, R.scaffold == bigchroms[15])
}

if (ref == "heradem") {
  spSummaryFilterSub = subset(spSummaryFilter, R.scaffold == "Herato1505")
}

ggplot(spSummaryFilterSub, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  # scale_x_continuous(limits = c(0,20)) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species")



## summarize all information ==================================================
buffer=25000
AllSummary = data.frame()
for (scaff in unique(spSummaryFilter$R.scaffold)) {
  s1 = subset(spSummaryFilter, R.scaffold == scaff)
  # s1 = subset(spSummaryFilter, R.scaffold == "Herato2101")
  s1 = s1[order(s1$Bsta),]
  s1$species = as.character(s1$species)
  s1$count = 1
  s1$spCnt = 1
  s1$eraCnt = 0
  s1$melCnt = 0
  s2$midCnt = 0
  s2 = s1[1,]
  for (rw in 2:nrow(s1)) {
    if (s1$Bsta[rw]-s2$Bsta[nrow(s2)] < buffer & s1$Bend[rw]-s2$Bend[nrow(s2)] < buffer) {
      s2$Bsta[nrow(s2)] = min(s1$Bsta[rw],s2$Bsta[nrow(s2)])
      s2$Bend[nrow(s2)] = max(s1$Bend[rw],s2$Bend[nrow(s2)])
      s2$count[nrow(s2)] = s2$count[nrow(s2)] + 1
      s2$spCnt[nrow(s2)] = s2$spCnt[nrow(s2)] + 1
      #
      if (s2$species[nrow(s2)] != s1$species[rw]) {
        s2$species[nrow(s2)] = paste0(s2$species[nrow(s2)],",",s1$species[rw])
      }
    }
    if (s1$Bsta[rw]-s2$Bsta[nrow(s2)] >= buffer | s1$Bend[rw]-s2$Bend[nrow(s2)] >= buffer) {
      s2 = rbind(s2, s1[rw,])
    }
  }
  for (i in 1:nrow(s2)) {
    #
    s2$species[i] = paste(unique(strsplit(s2$species[i], ",")[[1]]), collapse = ",")
    s2$spCnt[i] = length(unique(strsplit(s2$species[i], ",")[[1]]))
    # count species per clade
    sppList = unique(strsplit(s2$species[i], ",")[[1]]) 
    s2$melCnt[i] = length(sppList[sppList %in% mel_clade])
    s2$eraCnt[i] = length(sppList[sppList %in% era_clade])
    s2$midCnt[i] = length(sppList[sppList %in% mid_clade])
  }
  AllSummary = rbind(AllSummary, s2)
}

# filter
if (ref == "hmelv25") { AllSummaryFilter = subset(AllSummary, melCnt >= 2) }
if (ref == "heradem") { AllSummaryFilter = subset(AllSummary, eraCnt >= 2) }
# AllSummaryFilter = subset(AllSummary, spCnt >= 2)
filename=paste0("/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/vOrigin/",ref,".InvCandidates.2spp.txt")
write.table(file=filename, AllSummaryFilter, row.names = F, col.names=F, quote=F, sep="\t")

## cross datasets
dfFinal = data.frame()
for (rw in 1:nrow(AllSummaryFilter)) {
  ch = AllSummaryFilter$R.scaffold[rw]
  st = AllSummaryFilter$Bsta[rw]
  en = AllSummaryFilter$Bend[rw]
  s1 = subset(spSummaryFilter, R.scaffold == ch & Bsta >= st & Bend <= en)
  dfFinal = rbind(dfFinal, s1)
}
ggplot(dfFinal, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~R.scaffold, nrow=4)


# plot
ggplot(AllSummaryFilter, aes(x = (Bsta+Bend)/2000000, y=1, width=(Bend-Bsta)/1000000, height=0.8)) +
  geom_tile(col="black", alpha=0.5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~R.scaffold) + xlab("Chromosome Position(Mb)") + ylab("Species")

