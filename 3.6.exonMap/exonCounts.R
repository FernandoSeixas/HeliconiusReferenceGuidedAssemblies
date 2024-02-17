## load libraries  ============================================
require(tidyr)
require(ggplot2)
require(ggrepel)

## read data ==================================================
spp = "hhec"
chr = "hmel08"
cnb = "0008"
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/"
fil = paste0(chr, ".exons.2.",spp,".mm2")
inp = read.table(paste0(dir,fil), fill = TRUE)
names(inp) = c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
inp$propMap = inp$M.Bases/inp$Q.length
#plot(density(inp$propMap))
# separate first column into gene, exon and original coordinates
inpComplex = inp %>% separate(Q.scaffold, c("gene", "exonNb","chrom","sta","end"), sep="_")
# geneExon
inpComplex$geneExon = paste0(inpComplex$gene,"-",inpComplex$exonNb)
# filter by chromosome
inpFilter = subset(inpComplex, R.scaffold == paste0(spp,cnb))
subset(inpComplex, gene == "HMEL037926g1.t1")


## Exon Counts ==================================================
geneExonNames = read.table(paste0(dir,chr,".target.geneExons.txt"))[,1]
df = data.frame(geneExon = geneExonNames, all = rep(0, length(geneExonNames)), flt = rep(0, length(geneExonNames)))
head(df)
for (rw in 1:nrow(df)) {
  exon = df$geneExon[rw]
  df$all[rw] = nrow(subset(inpComplex, geneExon == exon))
  df$flt[rw] = nrow(subset(inpFilter, geneExon == exon))
}
df = df %>% separate(geneExon, c("gene", "exon"), sep="-")
df$geneExon = paste0(df$gene,"-",df$exon)
df$geneExon = factor(df$geneExon, levels = df$geneExon)
df$gene = factor(df$gene, levels = unique(df$gene))
# plot
ggplot(df) +
  geom_bar(aes(x=geneExon, y=-all), stat = "identity", fill="grey", col="darkgrey") +
  geom_bar(aes(x=geneExon, y=flt, fill=gene, col=gene), stat = "identity", alpha=0.25) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  geom_hline(yintercept = 1, alpha=0.5) +
  ggtitle(paste0(spp, "-", chr)) + ylim(-10,10)



# ## plot Exons along Chromosome ==================================================
# inpOrder = inpFilter[order(inpFilter$R.start),]
# inpOrder$ypos = 0
# step = 10
# inpOrder$ypos[seq(1,nrow(inpOrder),step)] = 0.1
# inpOrder$ypos[seq(2,nrow(inpOrder),step)] = 0.2
# inpOrder$ypos[seq(3,nrow(inpOrder),step)] = 0.3
# inpOrder$ypos[seq(4,nrow(inpOrder),step)] = 0.4
# inpOrder$ypos[seq(5,nrow(inpOrder),step)] = 0.5
# inpOrder$ypos[seq(6,nrow(inpOrder),step)] = 0.6
# inpOrder$ypos[seq(7,nrow(inpOrder),step)] = 0.7
# inpOrder$ypos[seq(8,nrow(inpOrder),step)] = 0.8
# inpOrder$ypos[seq(9,nrow(inpOrder),step)] = 0.9
# inpOrder$ypos[seq(10,nrow(inpOrder),step)] = 1.0
# # plot
# stp = 10000000
# sta = 4500000
# end = sta + stp
# ggplot(inpOrder, aes(x=(R.start+R.end)/2, y=ypos, label=exonNb)) +
#   geom_rect(aes(xmin=R.start-1000, xmax=R.end+1000, ymin=0, ymax=1, fill = gene)) +
#   xlim(min(inpOrder$R.start),9000000) + geom_text()
# 
# subset(inpOrder, R.start > 9000000)
# 
