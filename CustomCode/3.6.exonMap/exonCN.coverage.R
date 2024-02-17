## load libraries
require(stringr)
require(ggplot2)
require(data.table)
require(tidyr)

options(scipen = 999)

## read data
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/2.Cov2CN/1.ExonCov/"
lfiles = list.files(dir)
CovTable = data.frame()
for (fil in lfiles) {
  inp = read.table(paste0(dir,fil))
  names(inp) = c("scaffold","sta","end","readCount","coverage","ind")
  inp$geneExon = str_remove(fil,".coverage.txt")
  CovTable = rbind(CovTable, inp)
}
CovTable$geneExon = str_replace(CovTable$geneExon, "-RA", ".RA")
CovTable <- CovTable %>% separate(geneExon, c("gene","exon"), "-")
CovTable$geneExon = paste0(CovTable$gene,"-",CovTable$exon)

CovTable[grep(pattern = "RA", CovTable$geneExon),]


### Calculate relative coverage ==================================================
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/3.1.compareMap/1.2.rmdup/hmelv25/"
lfiles = list.files(dir, pattern = ".cov.w25s25.txt")
WinRefCoverage = data.frame()
for (fl in lfiles) {
  inp = read.table(paste0(dir,fl))
  names(inp) = c("scaffold","sta","end","readCount","coverage","ind")
  WinRefCoverage = rbind(WinRefCoverage, inp)
}
MedianRefCoverage = aggregate(WinRefCoverage$coverage, list(WinRefCoverage$ind), median)
names(MedianRefCoverage) = c("Individual","MedianCoverage")

# relative coverage
CovTable$relCov = 0
for (rw in 1:nrow(CovTable)) {
  ii = CovTable$ind[rw]
  CovTable$relCov[rw] = CovTable$coverage[rw] / MedianRefCoverage[MedianRefCoverage$Individual == ii,]$MedianCoverage
}

# transform species names
CovTable$ind = tolower(CovTable$ind)
CovTable$ind = str_replace(CovTable$ind, "hel_", "h")
CovTable$ind = str_replace(CovTable$ind, "lap_", "h")

### Plot =========================================================================
chr="hmel02"; scaff="Hmel202001o"
chr="hmel04"; scaff="Hmel204001o"
chr="hmel08"; scaff="Hmel208001o"
chr="hmel09"; scaff="Hmel209001o"
chr="hmel20"; scaff="Hmel220003o"
chr="hmel21"; scaff="Hmel221001o"
#
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/"
exons = read.table(paste0(dir,chr,".target.exons.txt"))[,1:2]
exons$V1 = as.character(exons$V1)
for (i in grep("RA", exons$V2)) {
  exons$V1[i] = as.character(paste0(exons$V1[i],".RA"))
  exons$V2[i] = str_remove(exons$V2[i], pattern = "RA-")
}
#
geneExons = paste0(exons$V1,"-",exons$V2)
subspp = c("hmel","hhec","hele","hpar")
xxx = subset(CovTable, scaffold == scaff & ind %in% subspp)
xxx$ind = factor(xxx$ind, levels = subspp)
xxx$gene = factor(xxx$gene, levels = unique(exons$V1) )
xxx$geneExon = factor(xxx$geneExon, levels = geneExons)
xxx = xxx[!is.na(xxx$geneExon),] # remove NA rows
# plot
ggplot(xxx) + 
  geom_bar(aes(x=geneExon, y=relCov, colour=gene, fill=gene), stat = "identity", alpha=1.0) +
  geom_hline(yintercept = 1, lty=2) +
  facet_wrap(~ind, nrow=4) + 
  theme(axis.text.x = element_blank()) +
  xlab("Exon") + ylab("Relative Coverage") +
  labs(fill = "Gene Transcript", colour = "Gene Transcript")

ggsave(paste0(dir,chr,"exonCN.eps"), dpi = 150, width = 313, height = 144, units = "mm")

