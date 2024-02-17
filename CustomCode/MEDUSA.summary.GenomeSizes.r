## load modules
require(stringr)

## variables
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/finalAssemblies/"
mel="hmelv25.chromLength.txt"
era="heradem.chromLength.txt"


## read data ==================================================
hmel=read.table(paste0(dir,mel))
hera=read.table(paste0(dir,era))
names(hmel) = c("scaffold", "length")
names(hera) = c("scaffold", "length")
# add species names
hmel$species = substr(hmel$scaffold,1,nchar(as.character(hmel$scaffold))-6)
hera$species = substr(hera$scaffold,1,nchar(as.character(hera$scaffold))-6)
# add chrom names
hmel$chrom = substr(hmel$scaffold,nchar(as.character(hmel$scaffold))-5, nchar(as.character(hmel$scaffold))-4)
hera$chrom = substr(hera$scaffold,nchar(as.character(hera$scaffold))-5, nchar(as.character(hera$scaffold))-4)


## sum same chromosome scaffold lengths =======================
sppList = read.table(paste0(dir,"sppList"))[,1]
chromTable = data.frame()
for (i in 1:length(sppList)) {
  spp = as.character(sppList[i])
  subMel = subset(hmel, species == spp)
  subEra = subset(hera, species == spp)
  dfMel = aggregate(subMel$length, list(subMel$chrom), sum); names(dfMel)=c("chrom","hmel.Size")
  dfEra = aggregate(subEra$length, list(subEra$chrom), sum); names(dfEra)=c("chrom","hera.Size")
  dfCmb = merge(dfMel, dfEra, by = "chrom")
  dfCmb$species = spp
  chromTable = rbind(chromTable, dfCmb)
}

melChromTable = dcast(chromTable[,c(4,1,2)], chrom ~ species)
eraChromTable = dcast(chromTable[,c(4,1,3)], chrom ~ species)



require(ggplot2)
require(ggrepel)
## compare CHROMOSOME sizes based on different guiding genomes
# correlation
ChromCorTest = data.frame()
for (spp in unique(chromTable$species)) {
  s1 = subset(chromTable, species == spp)
  ct = cor.test(s1$hmel.Size, s1$hera.Size, method = "spearman")
  df = data.frame(species = spp, rho = ct$estimate, pvalue = ct$p.value)
  ChromCorTest = rbind(ChromCorTest, df)
}
min(ChromCorTest$rho)
max(ChromCorTest$rho)
# plot
ggplot(chromTable) + 
  geom_point(aes(x=hmel.Size/1000000, y=hera.Size/1000000, colour=species), alpha=0.5) + 
  geom_abline(intercept = 0, slope = 1) +
  xlim(8,30) + ylim(8,30) +
  facet_wrap(~species) +
  xlab("Chromosome Size (Mb) - hmelv25 guided") +
  ylab("Chromosome Size (Mb) - heradem guided") +
  theme(aspect.ratio=1)

max(chromTable$hmel.Size)/1000000
max(chromTable$hera.Size)/1000000

## compare GENOME sizes based on different guiding genomes
hm = aggregate(chromTable$hmel.Size, list(chromTable$species), sum); names(hm) = c("species", "hmel.genomeSize")
he = aggregate(chromTable$hera.Size, list(chromTable$species), sum); names(he) = c("species", "hera.genomeSize")
cb = merge(hm, he, by="species")
# correlation test
ct = cor.test(cb$hmel.genomeSize, cb$hera.genomeSize, method = "spearman")
# plot
ggplot(cb, aes(x=hmel.genomeSize/1000000, hera.genomeSize/1000000, label=species)) + 
  geom_point() + geom_text_repel() +
  xlim(250,425) + ylim(250,425) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Genome Size (Mb) - hmelv25 guided") +
  ylab("Genome Size (Mb) - heradem guided") +
  theme(aspect.ratio=1)


## Genome Sizes by Group
# define groups
melClade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
midClade = c("hdor","hbur")
eraClade = c("hera","hhim","hhimfat","hsia","htel","hdem","hsar")
# assign to clade
cb$Clade = NA
cb$Clade = ifelse(cb$species %in% melClade, "mel", cb$Clade)
cb$Clade = ifelse(cb$species %in% midClade, "mid", cb$Clade)
cb$Clade = ifelse(cb$species %in% eraClade, "era", cb$Clade)
cb$Clade = factor(cb$Clade, levels = c("mel","mid","era"))
# # "BEST" genome size
# cb$GenomeSize = 0
# cb$GenomeSize = ifelse(cb$species %in% melClade, cb$hmel.genomeSize, cb$GenomeSize)
# plot
ggplot(cb) + 
  geom_boxplot(aes(x=Clade, y=hmel.genomeSize/1000000, fill=Clade), alpha=0.6) +
  geom_jitter(aes(x=Clade, y=hmel.genomeSize/1000000, color=Clade), alpha=0.6, pch=4) +
  scale_fill_manual(values = c("blue","orange","red")) +
  scale_color_manual(values = c("blue","orange","red")) +
  xlab("Clade") + ylab("Genome Size (Mb)")
ggplot(cb) + 
  geom_boxplot(aes(x=Clade, y=hera.genomeSize/1000000, fill=Clade), alpha=0.6) +
  geom_jitter(aes(x=Clade, y=hera.genomeSize/1000000, color=Clade), alpha=0.6, pch=4) +
  scale_fill_manual(values = c("blue","orange","red")) +
  scale_color_manual(values = c("blue","orange","red")) +
  xlab("Clade") + ylab("Genome Size (Mb)")

