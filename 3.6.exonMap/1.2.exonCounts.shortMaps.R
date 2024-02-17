## load libraries  ============================================
require(tidyr)
require(data.table)
require(ggplot2)
require(ggrepel)
require(stringr)

## hmel exon - short mappings to MEDUSA assemblies ============
# variables
chr = "hmel09"
cnb = str_replace(chr, "hmel","")
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.3.GenomeExpansion/exonMapping/1.Map2CN/0.minimap2/"
#lfiles = list.files(dir, pattern = ".mm2")
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
  # determine chromosome category
  chroms = unique(inpComplex$MainChrom)
  bigchr=c(); 
  for (nb in seq(1,9,1)) {nm = paste0(spp,"0",nb); bigchr = c(bigchr,nm) }
  for (nb in seq(10,21,1)) {nm = paste0(spp,nb); bigchr = c(bigchr,nm) }
  inpComplex$shortNameChrom = substr(inpComplex$MainChrom, 1, 6)
  inpComplex$ChromClass = NA
  inpComplex$ChromClass = ifelse(inpComplex$shortNameChrom %in% bigchr, "bigchr","nonchr")
  inpComplex$ChromClass = ifelse(inpComplex$shortNameChrom %in% paste0(spp, cnb), "target", inpComplex$ChromClass)
  # add to common data.frame
  dfExonMaps = rbind(dfExonMaps, inpComplex)
}

# filter mapping with adequate proportion of exon mapped
dfExonMapsFilter = subset(dfExonMaps, propMap >= 0.50)
# ggplot(dfExonMaps, aes(x=geneExon, y=propMap, colour=geneExon, fill=geneExon)) + 
#   geom_boxplot(alpha=0.1) +
#   geom_jitter(alpha=0.5) + 
#   facet_wrap(~spp) + ylim(0,1.2) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90))
# nrow(subset(dfExonMapsFilter, propMap >= 0.50))/nrow(dfExonMaps)



## Extend coordinates to then extract sequences ==============================
# extend coordinates
dfExonMapsFilter$extSta = NA
dfExonMapsFilter$extEnd = NA
dfExonMapsFilter$extSta = ifelse(dfExonMapsFilter$strand == "+", dfExonMapsFilter$R.start-(dfExonMapsFilter$Q.start)+1, dfExonMapsFilter$extSta)
dfExonMapsFilter$extEnd = ifelse(dfExonMapsFilter$strand == "+", dfExonMapsFilter$R.end+(dfExonMapsFilter$Q.length-dfExonMapsFilter$Q.end)-0, dfExonMapsFilter$extEnd)
dfExonMapsFilter$extSta = ifelse(dfExonMapsFilter$strand == "-", dfExonMapsFilter$R.start-(dfExonMapsFilter$Q.length-dfExonMapsFilter$Q.end)+1, dfExonMapsFilter$extSta)
dfExonMapsFilter$extEnd = ifelse(dfExonMapsFilter$strand == "-", dfExonMapsFilter$R.end+(dfExonMapsFilter$Q.start)+0, dfExonMapsFilter$extEnd)
dfExonMapsFilter$extSta = ifelse(dfExonMapsFilter$extSta < 0, 0, dfExonMapsFilter$extSta)
dfExonMapsFilter$extpropMap = abs(dfExonMapsFilter$extEnd-dfExonMapsFilter$extSta+1)/dfExonMapsFilter$Q.length
# write to file
write.table(file=paste0(dir,"../",chr,".mappedExons.sm.extended.txt"), dfExonMapsFilter, row.names=F, quote=F, sep="\t")
# subset to species of interest
# dfExonMapsFilter = subset(dfExonMapsFilter, spp %in% c("hmel","htim","hhec","hele","hpar"))


unique(dfExonMapsFilter$spp)

## Exon Counts ==================================================
geneExonNames = read.table(paste0(dir,"../../",chr,".target.geneExons.txt"))[,1]
geneExonNames = str_replace(geneExonNames, "\\-RA", "\\.RA")
dfCNV = data.frame()
sppList = c("hmel","htim","hcyd", "hbes", "hnum", "hhec","hele","hpar", "hera")
for (ex in geneExonNames) {
  for (sp in sppList) {
    st = subset(dfExonMapsFilter, spp == sp & geneExon == ex & ChromClass == "target")
    sb = subset(dfExonMapsFilter, spp == sp & geneExon == ex & ChromClass == "bigchr")
    sn = subset(dfExonMapsFilter, spp == sp & geneExon == ex & ChromClass == "nonchr")
    df = data.frame(geneExon = ex, tgt = nrow(st), big = nrow(sb)+nrow(st), non=nrow(sn)+nrow(sb)+nrow(st), spp=sp)
    dfCNV = rbind(dfCNV, df)
  }
}
dfCNV$geneExon = str_replace(dfCNV$geneExon, "\\-RA", "\\.RA")
dfCNV = dfCNV %>% separate(geneExon, c("gene", "exon"), sep="-")
dfCNV$geneExon = paste0(dfCNV$gene,"-",dfCNV$exon)
dfCNV$geneExon = factor(dfCNV$geneExon, levels = unique(dfCNV$geneExon))
dfCNV$gene = factor(dfCNV$gene, levels = unique(dfCNV$gene))
dfCNV$spp = factor(dfCNV$spp, levels=sppList)
# plot
dfCNVsub = subset(dfCNV, spp %in% c("hmel","hhec","hele","hpar"))
names(dfCNVsub)[1] = "GeneTranscript"

ggplot(dfCNVsub) +
  geom_bar(aes(x=geneExon, y=non), fill="grey", col="black", stat = "identity") +
  geom_bar(aes(x=geneExon, y=big), fill="black", col="black", stat = "identity") +
  geom_bar(aes(x=geneExon, y=tgt, fill=GeneTranscript, col=GeneTranscript), stat = "identity") +
  geom_hline(yintercept = 1, lty=2) +
  facet_wrap(~spp, nrow=5) +
  theme(axis.text.x = element_blank()) +
  # theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) + 
  # ggtitle(chr) + 
  xlab("Gene Transcript") + ylab("Number of Alignments")

ggsave(paste0(dir,chr,"exonMap.eps"), dpi = 150, width = 313, height = 144, units = "mm")


# mmm = subset(dfCNV, spp == "hhec")
# mmm = subset(dfCNV, spp == "hele")
# mmm = subset(dfCNV, spp == "hpar")
# aggregate(mmm$tgt, list(mmm$gene), max)
# 
# #
# dfExonMapsFilter$geneExon = factor(dfExonMapsFilter$geneExon, levels = unique(dfExonMapsFilter$geneExon))
# dfExonMapsFilter$gene = factor(dfExonMapsFilter$gene, levels = unique(dfExonMapsFilter$gene))
# dfExonMapsFilter$spp = factor(dfExonMapsFilter$spp, levels=c("hmel", "htim", "hhec", "hele", "hpar"))
# dfExonMapsFilter$ChromClass = factor(dfExonMapsFilter$ChromClass, levels = c("nonchr","bigchr","target"))
# dfExonMapsFilter$alpha = NA
# dfExonMapsFilter$alpha = ifelse(dfExonMapsFilter$ChromClass == "target", 0.30, dfExonMapsFilter$alpha)
# dfExonMapsFilter$alpha = ifelse(dfExonMapsFilter$ChromClass == "bigchr", 0.10, dfExonMapsFilter$alpha)
# dfExonMapsFilter$alpha = ifelse(dfExonMapsFilter$ChromClass == "nonchr", 0.20, dfExonMapsFilter$alpha)
# ggplot(dfExonMapsFilter) + 
#   geom_bar(aes(geneExon, fill=gene, col=ChromClass, alpha=alpha)) +
#   facet_wrap(~spp, nrow=5) +
#   scale_color_manual(values = c("black","black","black")) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) + 
#   xlab("Gene-Exon") + ylab("Copy Number (CN)") + ggtitle(chr)
# 
# ## CNV vs exon length
# dfCNV = data.frame()
# sppList = c("hhec","hele","hpar")
# for (ex in geneExonNames) {
#   for (sp in sppList) {
#     # sf = subset(dfExonMaps, spp == sp & geneExon == ex)
#     sf = subset(dfExonMapsFilter, spp == sp & geneExon == ex)
#     if (nrow(sf) == 0) {len=NA}
#     if (nrow(sf) > 0) {len=unique(sf$end)-unique(sf$sta)+1}
#     df = data.frame(geneExon = ex, CN = nrow(sf), spp=sp, len=len)
#     dfCNV = rbind(dfCNV, df)
#   }
# }
# dfCNV$lenbin = as.factor(round(((dfCNV$len/3)-floor(dfCNV$len/3))*3))
# dfCNV = dfCNV[!is.na(dfCNV$lenbin),]
# ggplot(dfCNV, aes(x=lenbin, y=CN, fill=spp, colour=spp)) + geom_boxplot(alpha=0.25) + geom_jitter(alpha=0.5) + facet_wrap(~spp)
# 
# dfCNV = dfCNV %>% separate(geneExon, c("gene","exon"),"-")
# head(dfCNV)
# 
# ## ======================================================================
# focalSpp = subset(dfExonMapsFilter, spp %in% c("hhec","hele","hpar"))
# focalSpp = focalSpp %>% separate(geneExon, c("gene","exon"), sep="-")
# focalSpp = focalSpp %>% separate(gene, c("gene","transcript"), sep="[.]")
# focalSpp = subset(focalSpp, transcript == "t1")
# focalSpp$geneExon = paste0(focalSpp$gene,".",focalSpp$transcript,"-",focalSpp$exon)
# 
# 
# ## gene Counts per Scaffold =========================
# geneCountsPerScaffold = data.frame()
# for (scaff in unique(focalSpp$NewScaff)) {
#   s1 = subset(focalSpp, NewScaff == scaff)
#   for (sp in unique(s1$spp)) {
#     s2 = subset(s1, spp == sp)
#     gCount = length(unique(s2$gene))
#     df1 = data.frame(spp = sp, scaffold = scaff, geneCount = gCount)
#     geneCountsPerScaffold = rbind(geneCountsPerScaffold, df1)
#   }
# }
# multiGenes = subset(geneCountsPerScaffold, geneCount > 1)
# table(multiGenes$spp)
# table(geneCountsPerScaffold$spp)
# 
# 
# ## gene complete/partial duplications ==============================
# mapsPerScaffold = data.frame()
# # analyze each SCAFFOLD separately
# for (scaff in unique(focalSpp$NewScaff)) {
#   s1 = subset(focalSpp, NewScaff == scaff) 
#   # analyze each SPECIES separately
#   for (sp in unique(s1$spp)) {
#     s2 = subset(s1, spp == sp)
#     ScaLen = unique(s2$R.length)
#     mGnCount = nrow(s2) # count how many times each gene occurs (may be different exons) 
#     mExCount = max(table(s2$geneExon))
#     df1 = data.frame(spp = sp, scaffold = scaff, scaffoldLen = ScaLen, mapCount = mGnCount, mexCount = mExCount)
#     mapsPerScaffold = rbind(mapsPerScaffold, df1)
#   }
# }
# table(subset(mapsPerScaffold, mapCount >= 2)$spp)
# table(subset(mapsPerScaffold, mexCount >= 2)$spp)
# mx = max(mapsPerScaffold$scaffoldLen)
# ggplot(mapsPerScaffold) + geom_histogram(aes(x=scaffoldLen/1000), breaks=seq(0,mx+10000,1000)/1000)
# 
# 
# ## focus on larger scaffolds
# bigScaffolds = subset(mapsPerScaffold, scaffoldLen >= 25000)
# table(subset(bigScaffolds, mapCount >= 2)$spp)
# table(subset(bigScaffolds, mexCount >= 2)$spp)
# hist(bigScaffolds$mapCount, breaks=seq(1,50,1))
# # plot candidate scaffold
# multiMaps = subset(bigScaffolds, mapCount >= 5)
# length(multiMaps$spp)
# i = 1
# # tgScaffold = subset(focalSpp, NewScaff == multiMaps$scaffold[i] & spp == multiMaps$spp[i])
# tgScaffold = subset(focalSpp, NewScaff == multiMaps$scaffold[i] & spp == multiMaps$spp[i] & Q.length >= 150)
# ggplot(tgScaffold) +
#   geom_rect(aes(xmin=R.start, xmax=R.end, ymin=0, ymax=1, fill=gene, col=strand), alpha=0.25) +
#   # geom_rect(aes(xmin=extSta, xmax=extEnd, ymin=0, ymax=1, fill=gene, col=strand), alpha=0.25) +
#   geom_text(aes(x=(R.start+R.end)/2, y=1.1, label=exonNb), angle=90) +
#   xlab("Scaffold Positon (bp)") +
#   # facet_wrap(~gene) +
#   ggtitle(paste0(multiMaps$spp[i],":",multiMaps$scaffold[i])) + xlab("Scaffold Position (bp)")
# 
# 
# 
# #################################################################################
# # define last exon per gene
# dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/"
# exonsTable = read.table(paste0(dir,chr,".target.exonsFrame.txt"), header=T)
# exonsTable = exonsTable %>% separate(gene, c("gene","transcript"),"[.]")
# exonsTable = subset(exonsTable, transcript == "t1")
# exonsTable$exonNb = str_remove(exonsTable$exon, "E") 
# # get last exon per gene
# lastExon = data.frame()
# for (gn in unique(exonsTable$gene)) {
#   s1=subset(exonsTable, gene == gn)
#   df = data.frame(gene = gn, exonNb = max(as.numeric(s1$exonNb)), lastExon = paste0("E",max(as.numeric(s1$exonNb))))
#   lastExon = rbind(lastExon, df)
# }
# lastExon = subset(lastExon, exonNb > 3) # subset to genes with more than 2 exons
# 
# ## focus on gene and chose scaffolds with the first and last exon
# candidates = data.frame()
# for (i in 1:nrow(lastExon)) {
#   focusGene = subset(focalSpp, gene == lastExon$gene[i])
#   exl = lastExon$lastExon[i]; exlNb = as.numeric(str_remove(exl, "E")); exp = paste0("E", (exlNb-1))
#   # analyze each SCAFFOLD separately
#   for (scaff in unique(focusGene$NewScaff)) {
#     s1 = subset(focusGene, NewScaff == scaff) 
#     # analyze each SPECIES separately
#     for (sp in unique(s1$spp)) {
#       s2 = subset(s1, spp == sp)
#       # fst = nrow(subset(s2, exon == "E1" | exon == "E2"))
#       # lst = nrow(subset(s2, exon == exl | exon == exp))
#       fst = nrow(subset(s2, exon == "E1"))
#       lst = nrow(subset(s2, exon == exl))
#       if (fst >= 1 & lst >= 1) { candidates = rbind(candidates, s2) }
#     }
#   }
# }
# ll = unique(paste0(candidates$spp,"-",candidates$NewScaff))
# candDF = data.frame()
# for (cand in ll) {
#   sp = strsplit(cand, "-")[[1]][1]
#   ns = strsplit(cand, "-")[[1]][2]
#   df = data.frame(spp = sp, NewScaff = ns)
#   candDF = rbind(candDF, df)
# }
# nrow(candDF)
# i = 7
# tgScaffold = subset(focalSpp, NewScaff == candDF$NewScaff[i] & spp == candDF$spp[i])
# sp=unique(tgScaffold$spp[i])
# ns=unique(tgScaffold$NewScaff[i])
# ch=unique(tgScaffold$MainChrom[i])
# ggplot(tgScaffold) +
#   geom_rect(aes(xmin=R.start, xmax=R.end, ymin=0, ymax=1, fill=gene, col=strand), alpha=0.25) +
#   # geom_rect(aes(xmin=extSta, xmax=extEnd, ymin=0, ymax=1, fill=gene, col=strand), alpha=0.25) +
#   geom_text(aes(x=(R.start+R.end)/2, y=1.1, label=exonNb), angle=90, nudge_y = -0.2) +
#   xlab("Scaffold Positon (bp)") +
#   facet_wrap(~gene, nrow=6) +
#   ggtitle(paste0(ch,">",ns,":",sp)) + xlab("Scaffold Position (bp)")
# 
# #################################################################################
# 
# 
# 
# # ## plot Exons along Chromosome ==================================================
# # dfExonMapsFilterOrder = dfExonMapsFilter[order(dfExonMapsFilter$R.start),]
# # dfExonMapsFilterOrder$ypos = 0
# # step = 10
# # dfExonMapsFilterOrder$ypos[seq(1,nrow(dfExonMapsFilterOrder),step)] = 0.1
# # dfExonMapsFilterOrder$ypos[seq(2,nrow(dfExonMapsFilterOrder),step)] = 0.2
# # dfExonMapsFilterOrder$ypos[seq(3,nrow(dfExonMapsFilterOrder),step)] = 0.3
# # dfExonMapsFilterOrder$ypos[seq(4,nrow(dfExonMapsFilterOrder),step)] = 0.4
# # dfExonMapsFilterOrder$ypos[seq(5,nrow(dfExonMapsFilterOrder),step)] = 0.5
# # dfExonMapsFilterOrder$ypos[seq(6,nrow(dfExonMapsFilterOrder),step)] = 0.6
# # dfExonMapsFilterOrder$ypos[seq(7,nrow(dfExonMapsFilterOrder),step)] = 0.7
# # dfExonMapsFilterOrder$ypos[seq(8,nrow(dfExonMapsFilterOrder),step)] = 0.8
# # dfExonMapsFilterOrder$ypos[seq(9,nrow(dfExonMapsFilterOrder),step)] = 0.9
# # dfExonMapsFilterOrder$ypos[seq(10,nrow(dfExonMapsFilterOrder),step)] = 1.0
# # # plot
# # stp = 300000
# # staP = 3100000
# # endP = staP + stp
# # ggplot(dfExonMapsFilterOrder, aes(x=(as.numeric(sta)+as.numeric(end))/2000000, y=ypos, label=exonNb)) +
# #   geom_rect(aes(xmin=(sta-1000)/1000000, xmax=(end+1000)/1000000, ymin=0, ymax=1, fill = gene)) +
# #   xlim(staP/1000000, endP/1000000) + geom_text() +
# #   facet_wrap(~spp)
