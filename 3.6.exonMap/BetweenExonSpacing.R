## load libraries  ============================================
require(tidyr)
require(data.table)
require(ggplot2)
require(ggrepel)
require(stringr)
require(ggrepel)

## read data ==================================================
# variables
chr = "hmel08"
cnb = str_replace(chr, "hmel", "00")
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/exon2indepScaffolds/"
lfiles = list.files(dir, pattern = ".lr.mm2")
lfiles = lfiles[lfiles %like% chr]
# read data
dfExonMaps = data.frame()
sppList = c()
for (fl in lfiles) {
  inp = read.table(paste0(dir,fl), fill = TRUE, stringsAsFactors = F)
  # names(inp) = c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2")
  names(inp) = c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
  inp$propMap = inp$Nb.Bases/inp$Q.length
  # separate first column into gene, exon and original coordinates
  inpComplex = inp %>% separate(Q.scaffold, c("gene", "exonNb","chrom","sta","end"), sep="_")
  inpComplex$sta = as.numeric(inpComplex$sta)
  inpComplex$end = as.numeric(inpComplex$end)
  # separate Ref chromosome names into original chrom/scaffold name and new scaffold name
  inpComplex = inpComplex %>% separate(R.scaffold, c("MainChrom", "NewScaff"), sep="-")
  # geneExon
  inpComplex$geneExon = paste0(inpComplex$gene,"-",inpComplex$exonNb)
  # add spp name
  # spp = str_remove(fl, pattern = ".splitScaffolds.mm2")
  spp = str_remove(fl, pattern = ".splitScaffolds.lr.mm2")
  spp = str_remove(spp, paste0(chr, ".exons.2."))
  inpComplex$spp = spp
  sppList = c(sppList, spp)
  # add to common data.frame
  dfExonMaps = rbind(dfExonMaps, inpComplex)
}
# filter for mappings to the expected chromosome
dfExonMapsFilter = subset(dfExonMaps, MainChrom %in% paste0(sppList, cnb))
# filter mapping with adequate proportion of exon mapped
# ggplot(dfExonMapsFilter, aes(x=geneExon, y=propMap, colour=geneExon, fill=geneExon)) + 
#   geom_boxplot(alpha=0.1) +
#   geom_jitter(alpha=0.5) + 
#   facet_wrap(~spp) + ylim(0,1.2) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90))
dfExonMapsFilter = subset(dfExonMapsFilter, propMap >= 0.50)
nrow(subset(dfExonMapsFilter, propMap >= 0.50))/nrow(dfExonMaps)

## Expected Exon Distances( based on melpomene reference genome)
# read exon coordinates file
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/"
exonCoords = read.table(paste0(dir,chr,".target.exons.txt"))
names(exonCoords) = c("gene","exon","scaffold","sta","end")
ExpExonDistances = data.frame()
for (gn in unique(exonCoords$gene)) {
  # chose exons in 1 gene
  s1 = subset(exonCoords, gene == gn)
  # if more than 1 exon then, for each exon in that gene calculate distance to other exons
  if (nrow(s1) > 1) {
    for (ex1 in 1:nrow(s1)) {
      en1 = s1$end[ex1]
      for (ex2 in ex1:nrow(s1)) {
        st2 = s1$sta[ex2]
        df = data.frame(gene = gn, exon.1 = s1$exon[ex1], exon.2 = s1$exon[ex2], dist = st2-en1)
        ExpExonDistances = rbind(ExpExonDistances, df)
      }
    }
  }
}
ExpExonDistances$spp = "hmelv25"
ExpExonDistances$absDist = abs(ExpExonDistances$dist)
ExpExonDistances$comp = paste0(ExpExonDistances$gene, ExpExonDistances$exon.1, ExpExonDistances$exon.2)
ExpExonDistances$ExonPair = paste0(ExpExonDistances$exon.1, ExpExonDistances$exon.2)


## Observed Exon Distances ==================================================
ObsExonDistances = data.frame()
for (spname in sppList) {
  print(spname)
  df1 = data.frame()
  obs = subset(dfExonMapsFilter, spp == spname)
  for (scaff in unique(obs$NewScaff)) {
    subScaff = subset(obs, NewScaff == scaff)
    # if more than one map
    if (nrow(subScaff) > 1) {
      for (gn in unique(subScaff$gene)) {
        subGene = subset(subScaff, gene == gn)
        # if more than 1 exon then, for each exon in that gene calculate distance to other exons
        if (nrow(subGene) > 1) {
          subGene = subGene[order(subGene$R.start),]
          for (ex in 2:nrow(subGene)) {
            en1 = as.numeric(subGene$R.end[(ex-1)])
            st2 = as.numeric(subGene$R.start[ex])
            df = data.frame(gene = gn, exon.1 = subGene$exon[(ex-1)], exon.2 = subGene$exon[ex], dist = st2-en1, scaff = unique(subGene$NewScaff))
            df1 = rbind(df1, df)
          }
        }
      }
    }
  }
  df1$spp = spname
  ObsExonDistances = rbind(ObsExonDistances, df1)
}

for (rw in 1:nrow(ObsExonDistances)) {
  exo1 = as.numeric(str_remove(ObsExonDistances$exon.1[rw], "E"))
  exo2 = as.numeric(str_remove(ObsExonDistances$exon.2[rw], "E"))
  if (exo1 > exo2) {
    ObsExonDistances$exon.1[rw] = paste0("E", exo2)
    ObsExonDistances$exon.2[rw] = paste0("E", exo1)
  }
}
# bad pairs
BadPairs = droplevels(subset(ObsExonDistances, exon.1 == exon.2))
# good pairs
ObsExonDistances = subset(ObsExonDistances, exon.1 != exon.2)
ObsExonDistances$absDist = abs(ObsExonDistances$dist)
ObsExonDistances$comp = paste0(ObsExonDistances$gene, ObsExonDistances$exon.1, ObsExonDistances$exon.2)
ObsExonDistances$ExonPair = paste0(ObsExonDistances$exon.1, ObsExonDistances$exon.2)

#####
ObsExonDistances$type="obs"
ExpExonDistances$type="exp"

validComp=c(); 
for (i in seq(2,30,1)) { validComp = c(validComp,paste0("E",(i-1),"E",i)) }

ObsExonDistances = subset(ObsExonDistances, ExonPair %in% validComp)
ExpExonDistances = subset(ExpExonDistances, ExonPair %in% validComp)

ObsExonDistances$expdist = 0
ObsExonDistances$reldist = 0
for (nr in 1:nrow(ObsExonDistances)) {
  cp = ObsExonDistances$comp[nr]
  expdist = ExpExonDistances[ExpExonDistances$comp == cp,]$absDist
  ObsExonDistances$expdist[nr] = expdist
  ObsExonDistances$reldist[nr] = ObsExonDistances$absDist[nr]/expdist
}

ObsExonDistances$gene = factor(ObsExonDistances$gene, levels=unique(ExpExonDistances$gene))
ObsExonDistances$spp = factor(ObsExonDistances$spp, levels = c("hmel","hbes","hhec","hele","hpar"))

# plot
ggplot(ObsExonDistances) +
  geom_point(aes(x=expdist, y=absDist, colour=gene), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) + ylim(0,30000) + xlim(0,30000) +
  facet_wrap(~spp) + ggtitle(chr) + theme(legend.position = "none") + 
  xlab("Expected Distance (bp)") + ylab("Observed Distance (bp)")

ggplot(ObsExonDistances) +
  geom_boxplot(aes(x=gene, y=reldist, fill=gene, col=gene), alpha=0.01) +
  geom_jitter(aes(x=gene, y=reldist, fill=gene, col=gene), alpha=0.5) +
  # geom_hline(yintercept = 1) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~spp, nrow=4) + ggtitle(chr) + ylab("Relative Exon Distance - Obs/Exp") + xlab("Gene")

ggplot(ObsExonDistances) +
  geom_boxplot(aes(x=gene, y=log10(reldist), fill=gene, col=gene), alpha=0.01) +
  geom_jitter(aes(x=gene, y=log10(reldist), fill=gene, col=gene), alpha=0.5) +
  # geom_hline(yintercept = 1) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~spp, nrow=4) + ggtitle(chr) + ylab("Relative Exon Distance - Log10(Obs/Exp)") + xlab("Gene")

ggplot(ObsExonDistances) +
  geom_boxplot(aes(x=comp, y=log10(reldist), fill=gene, col=gene), alpha=0.01) +
  geom_jitter(aes(x=comp, y=log10(reldist), fill=gene, col=gene), alpha=0.5) +
  # geom_hline(yintercept = 1) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~spp, nrow=4) + ggtitle(chr) + ylab("Relative Exon Distance - Log10(Obs/Exp)") + xlab("Exon Pair")
