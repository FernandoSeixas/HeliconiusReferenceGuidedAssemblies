## load libraries
require(stringr)
require(ggplot2)
require(data.table)
require(tidyr)
require(gridExtra)

options(scipen = 999)

## read data
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/3.Cov2CN.popdata/1.ExonCov/"
lfiles = list.files(dir)
CovTable = data.frame()
for (fil in lfiles) {
  inp = read.table(paste0(dir,fil))
  names(inp) = c("scaffold","sta","end","readCount","coverage","ind")
  inp$geneExon = str_remove(fil,".coverage.txt")
  CovTable = rbind(CovTable, inp)
}
# update gene/exon names
CovTable$geneExon = str_replace(CovTable$geneExon, "-RA", ".RA")
CovTable <- CovTable %>% separate(geneExon, c("gene","exon"), "-")
CovTable$geneExon = paste0(CovTable$gene,"-",CovTable$exon)
CovTable[grep(pattern = "RA", CovTable$geneExon),]
# add species population name
CovTable$spp = substr(CovTable$ind, 1, 4) 
CovTable$pop = substr(CovTable$ind, 1, 12)
CovTable$pop = ifelse(CovTable$pop ==  "Hhec.lei.pan", "Hhec.unk.pan", CovTable$pop)
CovTable$pop = ifelse(CovTable$pop ==  "Hhec.zul.pan", "Hhec.unk.pan", CovTable$pop)
CovTable$pop = ifelse(CovTable$pop ==  "Hpar.ore.ecu", "Hpar.nap.ecu", CovTable$pop)
poplevels = c("Hmel.vul.col", "Hhec.unk.pan","Hele.bar.gui","Hpar.ser.per",
              "Hmel.ama.per", "Hhec.fel.per","Hele.yur.per","Hpar.yur.per",
              "Hmel.agl.per", "Hele.nap.ecu","Hpar.nap.ecu")
CovTable$pop = factor(CovTable$pop, levels = poplevels)
CovTable$geo = ifelse(CovTable$pop %in% c("Hpar.yur.per", "Hele.yur.per", "Hhec.fel.per", "Hele.nap.ecu", "Hpar.nap.ecu", "Hmel.agl.per", "Hmel.ama.per"), "sym", "alo")

unique(CovTable$pop)



### Calculate relative coverage ==================================================
dir = "/n/holyscratch01/mallet_lab/fseixas/2.elepar/1.2.rmdup/wincov.rmdup/"
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
  ii = as.character(CovTable$ind[rw])
  CovTable$relCov[rw] = CovTable$coverage[rw] / MedianRefCoverage[MedianRefCoverage$Individual == ii,]$MedianCoverage
}

# coverage by Gene / aggregate max coverage [among exons] 
CovTableGene = aggregate(CovTable$relCov, list(paste0(CovTable$ind,"-",CovTable$gene,"-",CovTable$scaffold)), max)
names(CovTableGene) = c("indGene", "RelativeCoverage")
CovTableGene <- CovTableGene %>% separate(indGene, c("ind","gene","scaffold"), "-")
CovTableGene$pop = substr(CovTableGene$ind, 1,12)
CovTableGene$pop = ifelse(CovTableGene$pop ==  "Hhec.lei.pan", "Hhec.unk.pan", CovTableGene$pop)
CovTableGene$pop = ifelse(CovTableGene$pop ==  "Hhec.zul.pan", "Hhec.unk.pan", CovTableGene$pop)
CovTableGene$pop = ifelse(CovTableGene$pop ==  "Hpar.ore.ecu", "Hpar.nap.ecu", CovTableGene$pop)
CovTableGene$pop = factor(CovTableGene$pop, levels = poplevels)
CovTableGene$spp = substr(CovTableGene$pop, 1, 4)
CovTableGene$spp = factor(CovTableGene$spp, levels = c("Hmel", "Hhec", "Hele", "Hpar"))
CovTableGene$geo = ifelse(CovTableGene$pop %in% c("Hpar.yur.per", "Hele.yur.per", "Hhec.fel.per", "Hele.nap.ecu", "Hpar.nap.ecu", "Hmel.agl.per", "Hmel.ama.per"), "sym", "alo")



# gene against gene
g2g = subset(CovTableGene, gene %in% c("HMEL037926g1.t1","HMEL008474g1.t1")) # hif1a vs cadherin 
g2g = subset(CovTableGene, gene %in% c("HMEL007446g1.t1","HMEL008474g1.t1")) # Dhc64C vs cadherin ! min,med,max
g2g = subset(CovTableGene, gene %in% c("HMEL006651g1.t1","HMEL008474g1.t1")) # capicua vs cadherin ! min,med,max
g2g = subset(CovTableGene, gene %in% c("HMEL009427g1.t1","HMEL008474g1.t1")) # trabuco vs cadherin ! min,med,max
ggplot(g2g) +
  geom_line(aes(x=gene, y=RelativeCoverage, group=ind, col=geo)) +
  geom_hline(yintercept = 1, alpha=0.5, lty=2) +
  facet_wrap(~spp)

aaa = subset(CovTableGene, gene == "HMEL007446g1.t1")
bbb = subset(CovTableGene, gene == "HMEL037926g1.t1")
cor.test(aaa$RelativeCoverage, bbb$RelativeCoverage)

g2g = subset(CovTableGene, gene %in% c("HMEL009427g1.t1", "HMEL037926g1.t1", "HMEL007446g1.t1", "HMEL008474g1.t1")) # trabuco vs cadherin ! min,med,max
ggplot(g2g) + 
  geom_boxplot(aes(x=gene, y=RelativeCoverage, fill=geo, col=geo), alpha=0.6) +
  facet_wrap(~spp, nrow=2, scales = "free_y") +
  xlab("Gene") + ylab("Relative Coverage") + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))




### Plot =========================================================================
RepeatRegions = data.frame(
  chr = c("hmel02","hmel04","hmel08","hmel09"),
  sca = c("Hmel202001o","Hmel204001o","Hmel208001o","Hmel209001o"),
  code = c("A","B","C","D"),
  plotNumber = c(1,2,3,4)
)

hecColor = "#1b9e77"
eleColor = "#7570b3"
parColor = "#d95f02"
melColor = "black"

#
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/"
for (r in 1:4) {
  # define region to plot
  s1 = RepeatRegions[r,]
  chr = as.character(s1$chr)
  scaff = as.character(s1$sca)
  title = as.character(s1$code)
  pnumb = as.character(s1$plotNumber)
  #
  exons = read.table(paste0(dir,chr,".target.exons.txt"))[,1:2]
  exons$V1 = as.character(exons$V1)
  for (i in grep("RA", exons$V2)) {
    exons$V1[i] = as.character(paste0(exons$V1[i],".RA"))
    exons$V2[i] = str_remove(exons$V2[i], pattern = "RA-")
  }
  geneExons = paste0(exons$V1,"-",exons$V2)
  xxx = subset(CovTable, scaffold == scaff)
  xxx$gene = factor(xxx$gene, levels = unique(exons$V1) )
  xxx$geneExon = factor(xxx$geneExon, levels = geneExons)
  xxx = xxx[!is.na(xxx$gene),]

  ## plot by gene =========================
  mmm = aggregate(xxx$relCov, list(paste0(xxx$pop,"-",xxx$gene)), median)
  names(mmm) = c("popGene", "RelativeCoverage")
  mmm <- mmm %>% separate(popGene, c("pop","gene"), "-")
  mmm$gene = factor(mmm$gene, levels = unique(exons$V1) )
  mmm$pop = factor(mmm$pop, levels = poplevels)
  mmm$spp = substr(mmm$pop, 1, 4)
  mmm$spp = factor(mmm$spp, levels = c("Hhec", "Hele", "Hpar","Hmel"))
  mmm$geo = ifelse(mmm$pop %in% c("Hpar.yur.per", "Hele.yur.per", "Hhec.fel.per"), "sym1", "alo")
  mmm$geo = ifelse(mmm$pop %in% c("Hele.nap.ecu", "Hpar.nap.ecu"), "sym2", mmm$geo)
  mmm = mmm[mmm$spp != "Hmel",]
  p = ggplot(mmm) +
    geom_line(aes(x=gene, y=RelativeCoverage, group=pop, col=pop, lty=geo)) +
    geom_hline(yintercept = 1, alpha=1.0) +
    scale_color_manual(values = c(hecColor, eleColor, parColor, hecColor, eleColor, parColor, eleColor, parColor)) +
    facet_wrap(~spp, nrow=1) + ggtitle(chr) +
    theme(
      axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1), 
      plot.title = element_text(face="bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
      ) +
    xlab("Gene") + ylab("Relative Coverage") +
    ggtitle(title)
  assign(paste0("p", pnumb), p)
}
pf = grid.arrange(p1,p2,p3,p4, nrow=4)
ggsave(file=paste0(dir,"popCoverage.eps"), pf, dpi = "print", units = "mm", width = 210, height = 297)



# plot
# ggplot(xxx) + 
#   geom_boxplot(aes(x=geneExon, y=relCov, colour=gene, fill=gene), alpha=0.6) +
#   facet_wrap(~pop, nrow=2) + ggtitle(chr)

## plot by gene =========================
mmm = aggregate(xxx$relCov, list(paste0(xxx$pop,"-",xxx$gene)), median)
names(mmm) = c("popGene", "RelativeCoverage")
mmm <- mmm %>% separate(popGene, c("pop","gene"), "-")
mmm$gene = factor(mmm$gene, levels = unique(exons$V1) )
mmm$pop = factor(mmm$pop, levels = poplevels)
mmm$spp = substr(mmm$pop, 1, 4)
mmm$spp = factor(mmm$spp, levels = c("Hhec", "Hele", "Hpar","Hmel"))

mmm$geo = ifelse(mmm$pop %in% c("Hpar.yur.per", "Hele.yur.per", "Hhec.fel.per"), "sym1", "alo")
mmm$geo = ifelse(mmm$pop %in% c("Hele.nap.ecu", "Hpar.nap.ecu"), "sym2", mmm$geo)

hecColor = "#2ca25f"
eleColor = "#2b8cbe"
parColor = "#e34a33"
melColor = "black"
mmm = mmm[mmm$spp != "Hmel",]
ggplot(mmm) +
  # geom_point(aes(x=gene, y=RelativeCoverage, group=pop, col=pop), alpha=0.6) +
  geom_line(aes(x=gene, y=RelativeCoverage, group=pop, col=pop, lty=geo), alpha=0.6) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = c(hecColor, eleColor, parColor, hecColor, eleColor, parColor, eleColor, parColor)) +
  facet_wrap(~spp, nrow=1) + ggtitle(chr) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
  ggtitle(chr)
unique(mmm$pop)
# ggplot(mmm) +
#   geom_point(aes(x=gene, y=RelativeCoverage, colour=gene, fill=gene), alpha=0.6) +
#   geom_hline(yintercept = 1) +
#   facet_wrap(~pop, nrow=2) + ggtitle(chr) +
#   theme(axis.text.x = element_text(angle=90)) +
#   ggtitle(chr)


