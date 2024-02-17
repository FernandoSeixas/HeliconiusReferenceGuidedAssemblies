## load libraries
require(stringr)
require(ggplot2)
require(data.table)
require(tidyr)

## read data
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/3.Cov2CN-PopLevel/1.ExonCov/"
lfiles = list.files(dir)
CovTable = data.frame()
for (fil in lfiles) {
  inp = read.table(paste0(dir,fil))
  names(inp) = c("scaffold","sta","end","readCount","coverage","ind")
  inp$geneExon = str_remove(fil,".coverage.txt")
  CovTable = rbind(CovTable, inp)
}
CovTable <- CovTable %>% separate(geneExon, c("gene","exon"), "-")
CovTable$geneExon = paste0(CovTable$gene,"-",CovTable$exon)
## assign individuals to populations
inds = unique(CovTable$ind)
# define pops
CovTable$pop = str_replace(CovTable$ind, paste(c(".001",".002",".003",".004",".005"), collapse="|"), "")
# reorder exons
CovTable$exon = str_remove(CovTable$exon, "E")
CovTable$exon = as.numeric(CovTable$exon)
CovTable = CovTable[with(CovTable, order(gene, exon)),]
# remove bad individual
CovTable = subset(CovTable, ind != "Hpar.pmd.per.003")

# define groups
ele = sort(inds[grepl(paste(c("Hele"), collapse = "|"), inds)])
par = sort(inds[grepl(paste(c("Hpar"), collapse = "|"), inds)])
elePar = sort(inds[grepl(paste(c("Hele.bar","Hele.tum"), collapse = "|"), inds)])
parSer = sort(inds[grepl(paste(c("Hpar.ser"), collapse = "|"), inds)])
eleAma = ele[!(ele %in% elePar)]
parAma = par[!(par %in% parSer)]
melAma = sort(inds[grepl(paste(c("Hmel.ama"), collapse = "|"), inds)])
melAgl = sort(inds[grepl(paste(c("Hmel.agl"), collapse = "|"), inds)])
ethAer = sort(inds[grepl(paste(c("Heth"), collapse = "|"), inds)])
hecFel = sort(inds[grepl(paste(c("Hhec"), collapse = "|"), inds)])
out = c(as.character(melAma), as.character(melAgl), as.character(ethAer), as.character(hecFel))
#
CovTable$Group = NA
CovTable$Group = ifelse(CovTable$ind %in% elePar, "eleAlo", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% parSer, "parAlo", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% eleAma, "eleAma", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% parAma, "parAma", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% melAma, "melAma", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% melAgl, "melAgl", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% ethAer, "ethAer", CovTable$Group)
CovTable$Group = ifelse(CovTable$ind %in% hecFel, "hecFel", CovTable$Group)



### Calculate relative coverage ==================================================
dir = "/n/scratchlfs/mallet_lab/fseixas/2.elepar/1.1.sort/"
lfiles = list.files(dir, pattern = ".coverage.w25s25.txt")
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
ggplot(CovTable) + geom_density(aes(x=log10(relCov), fill=pop)) + 
  facet_wrap(~Group) + theme(legend.position = "none")

### Plot =========================================================================
chr="hmel08"
scaff="Hmel208001o"
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/3.Cov2CN-PopLevel/"
#
exons = read.table(paste0(dir,"../",chr,".target.exons.txt"))[,1:2]
geneExons = paste0(exons$V1,"-",exons$V2)
# xxx = subset(CovTable, Group != "melAma" & Group != "melAgl" & Group != "ethAer" & Group != "hecFel" & gene == "HMEL037926g1.t1")
xxx = subset(CovTable, Group != "melAma" & Group != "melAgl" & scaffold == scaff)
xxx$geneExon = factor(xxx$geneExon, levels = geneExons)
xxx$Group = factor(xxx$Group, levels = c("parSer","elePar","parAma","eleAma","hecFel","ethAer"))
ggplot(xxx) + 
  geom_boxplot(aes(x=geneExon, y=relCov, colour=gene), alpha=0.2) +
  geom_hline(yintercept = 1) +
  facet_wrap(~Group, nrow=5) + 
  theme(axis.text.x = element_text(angle=90), legend.position = "none")


## Gene/Gene correlations ==================================================
# IngTable = subset(CovTable, !(ind %in% out))
IngTable = subset(CovTable, Group %in% c("eleAlo","eleAma","parAma","parAlo","hecFel"))
gene1 = "HMEL006651g1.t1"
gene2 = "HMEL002509g1.t1"
# gene1 = "HMEL037926g1.t1"
# gene2 = "HMEL008474g1.t1"
g1 = subset(IngTable, gene == gene1)
g2 = subset(IngTable, gene == gene2)
mg1 = max(g1$relCov)
mg2 = max(g2$relCov)
ggplot(g1) + geom_line(aes(x=exon, y=relCov, colour=pop, group=ind)) + facet_wrap(~Group, ncol=5) + scale_y_continuous(limits=c(0,mg1), breaks = seq(0,500,5))
ggplot(g2) + geom_line(aes(x=exon, y=relCov, colour=pop, group=ind)) + facet_wrap(~Group, ncol=5) + scale_y_continuous(limits=c(0,mg2), breaks = seq(0,500,5))
## median pops
p1 = data.frame()
for (ex1 in unique(g1$exon)) {
  s1=subset(g1, exon == ex1)
  a1=aggregate(s1$relCov, list(s1$Group), median)
  names(a1) = c("Group", "CN")
  a1$exon = ex1
  p1 = rbind(p1, a1)
}
p2 = data.frame()
for (ex2 in unique(g2$exon)) {
  s1=subset(g2, exon == ex2)
  a1=aggregate(s1$relCov, list(s1$Group), median)
  names(a1) = c("Group", "CN")
  a1$exon = ex2
  p2 = rbind(p2, a1)
}
mg1 = max(p1$CN)
mg2 = max(p2$CN)
ggplot(p1) + geom_line(aes(x=exon, y=CN, colour=Group)) + geom_point(aes(x=exon, y=CN, colour=Group)) + scale_y_continuous(limits=c(0,mg1), breaks = seq(0,500,5)) + ggtitle(gene1)
ggplot(p2) + geom_line(aes(x=exon, y=CN, colour=Group)) + geom_point(aes(x=exon, y=CN, colour=Group)) + scale_y_continuous(limits=c(0,mg2), breaks = seq(0,500,5)) + ggtitle(gene2)


## 
RhoTable = data.frame(matrix(nrow=length(unique(g1$geneExon)), ncol=length(unique(g2$geneExon))))
PvaTable = data.frame(matrix(nrow=length(unique(g1$geneExon)), ncol=length(unique(g2$geneExon))))
rownames(RhoTable) = unique(g1$geneExon); colnames(RhoTable) = unique(g2$geneExon) 
rownames(PvaTable) = unique(g1$geneExon); colnames(PvaTable) = unique(g2$geneExon) 
for (ex1 in 1:length(unique(g1$geneExon))) {
  e1 = unique(g1$geneExon)[ex1]
  subE1 = subset(g1, geneExon == e1)
  smpE1 = subE1[,c(6,12)]
  for (ex2 in 1:length(unique(g2$geneExon))) {
    e2 = unique(g2$geneExon)[ex2]
    subE2 = subset(g2, geneExon == e2)
    smpE2 = subE2[,c(6,12)]
    mrg = merge(smpE1,smpE2, by = "ind")
    ct = cor.test(mrg$relCov.x, mrg$relCov.y, method="spearman")
    RhoTable[ex1,ex2] = ct$estimate
    PvaTable[ex1,ex2] = ct$p.value
  }
}
heatmap(as.matrix(RhoTable), scale = "none")

#
subG1 = subset(g1, exon == "1")[,c(6,10,11,12)]
subG2 = subset(g2, exon == "23")[,c(6,10,11,12)]
subMr = merge(subG1,subG2, by = c("ind","pop","Group"))
subMr = subset(subMr, !(ind %in% out))
names(subMr) = c("ind","pop","Group","CN1","CN2")
ggplot(subMr, aes(x=CN1, y=CN2)) + geom_point(aes(colour=pop), alpha=0.50) + geom_smooth(method = "lm") + facet_wrap(~Group)
subMr = subset(subMr, Group == "parAma")
cor.test(subMr$CN1, subMr$CN2, method="spearman")

