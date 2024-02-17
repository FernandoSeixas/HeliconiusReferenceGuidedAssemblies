## load libraries
require(ggplot2)
require(stringr)

## read annotations ==================================================
# genes
gdir="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmel2.5/"
gfil="hmelv25.genes_coord.tab"
gene=read.table(paste0(gdir,gfil))
names(gene) = c("gene","scaffold","sta","end","strand")
# exons
edir="/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.2.compareMap/"
efil="hmel25_exons.unique.txt"
exon=read.table(paste0(edir,efil))
names(exon) = c("gene","scaffold","sta","end")
# repeats
rdir="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmel2.5/"
rfil="hmelv25_repeats.txt"
reps=read.table(paste0(rdir,rfil))
names(reps) = c("scaffold","sta","end")


## read coverage ==================================================
dir="/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.2.compareMap/" 
reg="chr20.1M-2M"
hele = read.table(paste0(dir,"hele.",reg,".coverage"), header=T)
hele$relcov = hele$COV/33

## plot data  =====================================================
# subregion
ch = "Hmel220003o"
st = 1050000
en = 1060000
# subset to region
subgene = subset(gene, scaffold == ch & end >= st & sta <= en)
subexon = subset(exon, scaffold == ch & end >= st & sta <= en)
subreps = subset(reps, scaffold == ch & end >= st & sta <= en)
subcove = subset(hele, REF == ch & POS >= st & POS <= en)
# gene name positions
subgene$ypos = 0
subgene[seq(1,nrow(subgene),4),]$ypos = -1
subgene[seq(2,nrow(subgene),4),]$ypos = -2
subgene[seq(3,nrow(subgene),4),]$ypos = -3
subgene[seq(4,nrow(subgene),4),]$ypos = -4
# plot
ggplot() +
  geom_line(data=subcove, aes(x=POS/1000000, y=relcov)) + 
  #geom_rect(data=subgene, aes(xmin=sta/1000000, xmax=end/1000000, ymin=-3, ymax=-2), fill="red", alpha=0.8) +
  #geom_rect(data=subexon, aes(xmin=sta/1000000, xmax=end/1000000, ymin=-5, ymax=0), fill="red", alpha=0.8) +
  geom_rect(data=subreps, aes(xmin=sta/1000000, xmax=end/1000000, ymin=-10, ymax=-5), fill="blue", alpha=0.8) +
  #geom_text(data=subgene, aes(x=(sta+end)/2000000, y=-11+(ypos*5), label=gene)) +
  xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") +
  ggtitle(ch)

subreps
