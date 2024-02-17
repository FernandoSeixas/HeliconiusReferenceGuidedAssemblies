## load modules
require(ggplot2)
require(ggrepel)
require(stringr)

## read data ==================================================
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
hmelv25 = read.table(paste0(dir,"InvCandidates.BreakpointsFilter-2-hmelv25.txt"), header=T)
heradem = read.table(paste0(dir,"InvCandidates.BreakpointsFilter-2-heradem.txt"), header=T)
heralat = read.table(paste0(dir,"InvCandidates.BreakpointsFilter-2-heralat.txt"), header=T)

# inversion lengths
hmelv25$Blen=hmelv25$Bend-hmelv25$Bsta+1 
heradem$Blen=heradem$Bend-heradem$Bsta+1 
heralat$Blen=heralat$Bend-heralat$Bsta+1 
# inversions passing length filter
nrow(subset(hmelv25, Blen > 50000))/nrow(hmelv25)
nrow(subset(heradem, Blen > 50000))/nrow(heradem)
nrow(subset(heralat, Blen > 50000))/nrow(heralat)


ggplot(hmelv25) + geom_histogram(aes(Blen/1000)) + xlim(0,1000)
ggplot(heradem) + geom_histogram(aes(Blen/1000)) + xlim(0,1000)
ggplot(heralat) + geom_histogram(aes(Blen/1000)) + xlim(0,1000)


## intersection of inversions supported by mapping to more than one reference
hmelv25$sppScaff = paste0(hmelv25$species,"-",hmelv25$Q.scaffold)
heradem$sppScaff = paste0(heradem$species,"-",heradem$Q.scaffold)
heralat$sppScaff = paste0(heralat$species,"-",heralat$Q.scaffold)

require(venn)
v3 = venn(list(hmelv25=hmelv25$sppScaff, heradem=heradem$sppScaff, heralat=heralat$sppScaff))


require(tidyr)
i3 = data.frame(sppScaff = intersect(intersect(hmelv25$sppScaff, heradem$sppScaff),heralat$sppScaff))
i3Merge = merge(i3, hmelv25, by="sppScaff")
i3Merge <- i3Merge %>% separate(sppScaff, c("species","Q.scaffold"))
i3Merge <- i3Merge[with(i3Merge, order(R.scaffold, Bsta)),]
# subset(i3Merge, R.scaffold == unique(i3Merge$R.scaffold)[20])




# merge inversions supported by the same scaffold in the same species
cmb_hmelv25_heradem = merge(hmelv25, heradem, by = c("species", "Q.scaffold"))
cmb_hmelv25_heralat = merge(hmelv25, heralat, by = c("species", "Q.scaffold"))
cmb_heradem_heralat = merge(heradem, heralat, by = c("species", "Q.scaffold"))
nrow(cmb_hmelv25_heradem)
nrow(cmb_hmelv25_heralat)
nrow(cmb_heradem_heralat)


cmb$region = paste0(cmb$R.scaffold.x,".",cmb$Bsta.x,".",cmb$Bend.x)

subset(cmb, R.scaffold.x == "Hmel215003o" & Bsta.x > 1000000 & Bsta.x < 2000000)[,c(1,3,4)]

## filter by coverage =========================================
dir="/n/scratchlfs/mallet_lab/fseixas/ziheng_yang/0.data/genome.stats/winCoverage/"
filenames <- list.files(dir, pattern="*.tohmelv25.cov.w25s25.txt", full.names=TRUE)
covTable = data.frame()
for (file in filenames) { inp = read.table(file); covTable = rbind(covTable, inp) }
names(covTable) = c("scaffold","start","end","sites","coverage","individual")
# transform individual/species names
covTable$individual = tolower(covTable$individual)
covTable$individual = str_replace(covTable$individual, "_","")
covTable$individual = str_replace(covTable$individual, "hel","h")
covTable$individual = str_replace(covTable$individual, "heratog","hera")
covTable$individual = str_replace(covTable$individual, "hhim_fat","hhimfat")
covTable$individual = str_replace(covTable$individual, "hhec_old","hhecold")
covTable$individual = str_replace(covTable$individual, "lapdor","hdor")
# retain only relevant individuals
covTable = subset(covTable, individual %in% unique(mel$species))
# get relative coverage by individual
covTable$relCoverage = 0
for (sp in unique(covTable$individual)) {
  maxCov = max(subset(covTable, individual == sp)$coverage)
  covTable$relCoverage = ifelse(covTable$individual == sp, covTable$coverage/maxCov, covTable$relCoverage)
}
# coverage outliers
covMedian = aggregate(covTable$coverage, list(covTable$individual), median)
names(covMedian) = c("species", "medianCov")
# define coverage outliers
covTable$outColor = alpha("black", 0)
for (spp in unique(covMedian$species)) {
  maxCov = subset(covMedian, species == spp)$medianCov*1.5
  covTable$outColor = ifelse(covTable$individual == spp & covTable$coverage >= maxCov, alpha("blue",0.8), covTable$outColor)
}
# assign coverage to breakpoints
cmb$staRelCovHmel = 0
cmb$endRelCovHmel = 0
for (rw in 1:nrow(cmb)) {
  # parameters
  sp = as.character(cmb$species[rw])
  ch = as.character(cmb$R.scaffold.x[rw])
  st = cmb$Bsta.x[rw]
  en = cmb$Bend.x[rw]
  # get relative coverage
  cmb$staRelCovHmel[rw] = subset(covTable, individual == sp & scaffold == ch & start < st & end > st)$coverage / median(subset(covTable, individual == sp)$coverage)
  cmb$endRelCovHmel[rw] = subset(covTable, individual == sp & scaffold == ch & start < en & end > en)$coverage / median(subset(covTable, individual == sp)$coverage)
}
# filter inversions if excessive coverage in breakpoints
cmbCovFilter = subset(cmb, staRelCovHmel < 1.25 & endRelCovHmel < 1.25)
cmbCovFilter = cmbCovFilter[with(cmbCovFilter, order(R.scaffold.x,Bsta.x)),]
nrow(cmbCovFilter)
subset(cmbCovFilter, species == "hele")

subset(cmbCovFilter, R.scaffold.x == "Hmel221001o")

## plot Inversions ==================================================
chr = "Hmel202001o"
# assign to clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")
# subset inversions
# ppp=subset(cmbCovFilter, R.scaffold.x == chr)
ppp=subset(cmb, R.scaffold.x == chr)
ppp$species = factor(ppp$species, levels = c(mid_clade, mel_clade, era_clade))
xmax=max(ppp$Bend.x + 1)/1000000
ppp$ymin = 0
ppp$ymax = 1
# subset coverage
covSub = subset(covTable, scaffold == chr)
covSub$sppPos = 0
ppp$sppPos = 0
for (i in 1:length(levels(ppp$species))) {
  sp = levels(ppp$species)[i]
  covSub$sppPos = ifelse(covSub$individual == sp, i, covSub$sppPos)
  ppp$sppPos = ifelse(ppp$species == sp, i, ppp$sppPos)
}
labels = data.frame()
for (i in 1:length(levels(ppp$species))) {
  sp = levels(ppp$species)[i]
  df = data.frame(species=sp, sppPos = i)
  labels=rbind(labels, df)
}
# plot
ggplot() +
  geom_line(data=covSub, aes(x=(start+end)/2000000, y=sppPos+(relCoverage*5), group=individual)) +
  geom_point(data=covSub, aes(x=(start+end)/2000000, y=sppPos+(relCoverage*5), group=individual), col=covSub$outColor) +
  geom_text(data=labels, aes(x=0, y=sppPos, label=species), hjust = 1.1) +
  geom_tile(data=ppp,aes(x = (Bsta.x+Bend.x)/2000000, y=sppPos, width=(Bend.x-Bsta.x)/1000000, height=0.8, fill=species), col="black", alpha=0.5) +
  ggtitle(chr) + xlab("Chromosome Position(Mb)") + ylab("Species") + scale_x_continuous(breaks = seq(0,15,1))

subset(cmbCovFilter, R.scaffold.x == "Hmel211001o")


## PCA =======================================================
spp = "hpar"
sppCan = subset(cmbCovFilter, species == spp)
sppCanSort = sppCan[with(sppCan, order(R.scaffold.x,Bsta.x)),]
sppCanSort$Bend.x-sppCanSort$Bsta.x

region = "Hmel201001o.1689638.1860873"
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.3.structVar/"
vec = read.table(paste0(dir,"4.pcaCandidates/",region,".eigenvec"))
val = read.table(paste0(dir,"4.pcaCandidates/",region,".eigenval"))
vec=vec[,2:ncol(vec)]
names(vec) = c("IND", paste0("PC",seq(1:(ncol(vec)-1))))
# define pops
elePop = vec$IND[grep("ele", vec$IND)]
parPop = vec$IND[grep("par", vec$IND)]
vec$POP = ifelse(vec$IND %in% elePop, "ele", "par")
# plot
ggplot(vec, aes(x=PC1, y=PC2, color=POP, label=IND)) + geom_point() +
  scale_color_manual(values = c("blue","red")) +
  # geom_text_repel(alpha=0.3) +
  theme_bw() + theme(panel.grid = element_blank()) + ggtitle(paste0(spp,"-",region)) +
  xlab(paste0("PC1 (",round(val$V1[1],2),"%)")) +
  ylab(paste0("PC2 (",round(val$V1[2],2),"%)"))
# plot(x=seq(1:nrow(vec)), y=vec$PC1, col=ifelse(vec$POP == "ele", "blue","red"))



# plot each region
for (i in 1:nrow(sppCanSort)) {
  region = sppCanSort$region[nrow(sppCanSort)]
  vec = read.table(paste0(dir,"4.pcaCandidates/",region,".eigenvec"))
  val = read.table(paste0(dir,"4.pcaCandidates/",region,".eigenval"))
  vec=vec[,2:ncol(vec)]
  names(vec) = c("IND", paste0("PC",seq(1:(ncol(vec)-1))))
  # define pops
  elePop = vec$IND[grep("ele", vec$IND)]
  parPop = vec$IND[grep("par", vec$IND)]
  vec$POP = ifelse(vec$IND %in% elePop, "ele", "par")
  # plot
  ggplot(vec, aes(x=PC1, y=PC2, color=POP, label=IND)) + geom_point() +
    scale_color_manual(values = c("blue","red")) +
    geom_text_repel(alpha=0.3) +
    theme_bw() + theme(panel.grid = element_blank()) + ggtitle(paste0(spp,"-",region)) +
    xlab(paste0("PC1 (",round(val$V1[1],2),"%)")) +
    ylab(paste0("PC2 (",round(val$V1[2],2),"%)"))
}

plot(vec$PC2, col=ifelse(vec$POP == "ele", "blue", "red"), pch=ifelse(vec$POP == "ele", 2, 3), ylab=(paste0("PC1 (",round(val$V1[1],2),"%)")), main=paste0(spp,"-",region))
#ggplot(vec, aes(x=POP, y=PC1, color=POP, label=IND)) + geom_boxplot() + geom_jitter()

region="Hmel209001o.5247342.5274952"