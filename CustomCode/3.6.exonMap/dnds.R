require(seqinr)
require(stringr)
require(data.table)
require(tidyr)
require(ggplot2)

##
dir1 = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/exon2indepScaffolds/5.frmseq/"
dir2 = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/exon2indepScaffolds/6.popseq/"

## All Pops Dn/Ds  ==================================================
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/"
chrom = "hmel08"
exonsTable = read.table(paste0(dir,chrom,".target.exonsFrame.txt"), header=T)

dfAllPop = data.frame()
for (gn in unique(exonsTable$gene)) {
  s1 = subset(exonsTable, gene == gn)
  for (ex in s1$exon) {
    print(paste0(gn,ex))
    lfiles = list.files(dir1, pattern = gn)
    lfiles = lfiles[lfiles %like% paste0(".",ex,".frame")]
    if (length(lfiles) > 0) {
      # check if files have data
      all.info = file.info(paste0(dir1,lfiles))$size
      # read fasta files
      if (all.info > 0) {all = read.alignment(paste0(dir1,lfiles), "fasta", forceToLower = TRUE) }
      # calculate pairwise dnds within pops and get mean
      all.md.dnds = NaN
      allint = (length(str_split(all$seq[[1]], "")[[1]])/3) / (floor(length(str_split(all$seq[[1]], "")[[1]])/3))
      if (all$nb > 1 & allint == 1) {
        all.dnds = kaks(all, rmgap = T); 
        all.md.dn = mean(as.vector(all.dnds$ka), na.rm = T);
        all.md.ds = mean(as.vector(all.dnds$ks), na.rm = T);
        all.md.dnds = mean(as.vector(all.dnds$ka/all.dnds$ks), na.rm = T);
      }
      #
      df1 = data.frame(gene = gn, exon = ex, all.dn = all.md.dn, all.ds = all.md.ds, all.dnds = all.md.dnds)
      dfAllPop = rbind(dfAllPop, df1)
    }
  }
}
#
dfAllPop$geneExon = paste0(dfAllPop$gene,dfAllPop$exon)
dfAllPop$geneExon = factor(dfAllPop$geneExon, levels = unique(dfAllPop$geneExon))
dfAllPop$gene = factor(dfAllPop$gene, levels = unique(dfAllPop$gene))
# transform negative dn/ds values to 0 and infinite values to NA
dfAllPop$all.dnds2 = dfAllPop$all.dn/dfAllPop$all.ds
#dfAllPop$all.dnds = ifelse(dfAllPop$all.dnds < 0, 0, dfAllPop$all.dnds)
#dfAllPop$all.dnds[!is.finite(dfAllPop$all.dnds)] = NA
# plot dn/ds´
mx = ifelse(max(dfAllPop$all.dnds2, na.rm=T) > 2, max(dfAllPop$all.dnds2, na.rm=T), 2)
# ggplot(dfAllPop, aes(x=gene, y=all, colour=gene)) + 
#   geom_point(alpha=0.75) + ylim(0,mx) +
#   theme(axis.text.x = element_text(angle=90), legend.position = "none") 
ggplot(dfAllPop) +
  geom_point(aes(x=geneExon, y=-all.dn, fill=gene), col="black", alpha=0.50) + 
  geom_point(aes(x=geneExon, y=-all.ds, fill=gene), col="blue", alpha=0.50) + 
  geom_bar(aes(x=geneExon, y=all.dnds2, colour=gene,fill=gene), alpha=0.25, stat="identity") + 
  ylim(-0.25,mx) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5), legend.position = "none")



## Within Pops Dn/Ds ==================================================
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.6.exonMapping/"
chrom = "hmel08"
exonsTable = read.table(paste0(dir,chrom,".target.exonsFrame.txt"), header=T)

dfWithinPop = data.frame()
for (gn in unique(exonsTable$gene)) {
  s1 = subset(exonsTable, gene == gn)
  for (ex in s1$exon) {
    print(paste0(gn,".",ex))
    lfiles1 = list.files(dir1, pattern = gn); lfiles1 = lfiles1[lfiles1 %like% paste0(".",ex,".frame")]
    lfiles = list.files(dir2, pattern = gn); lfiles  = lfiles[lfiles %like% paste0(".",ex,".h")]
    if (length(lfiles) > 0) {
      # check if files have data
      hall.info = file.info(paste0(dir1,lfiles1[lfiles1 %like% ".frame.fasta"]))$size
      trio.info = file.info(paste0(dir2,lfiles[lfiles %like% ".htri.popalign.fasta"]))$size
      hhec.info = file.info(paste0(dir2,lfiles[lfiles %like% ".hhec.popalign.fasta"]))$size
      hele.info = file.info(paste0(dir2,lfiles[lfiles %like% ".hele.popalign.fasta"]))$size
      hpar.info = file.info(paste0(dir2,lfiles[lfiles %like% ".hpar.popalign.fasta"]))$size
      hmel.info = file.info(paste0(dir2,lfiles[lfiles %like% ".hmelOut.popalign.fasta"]))$size
      # read fasta files
      if (hall.info > 0) {hall = read.alignment(paste0(dir1,lfiles1[lfiles1 %like% ".frame.fasta"]), "fasta", forceToLower = TRUE) }
      if (trio.info > 0) {trio = read.alignment(paste0(dir2,lfiles[lfiles %like% ".htri.popalign.fasta"]), "fasta", forceToLower = TRUE) }
      if (hhec.info > 0) {hhec = read.alignment(paste0(dir2,lfiles[lfiles %like% ".hhec.popalign.fasta"]), "fasta", forceToLower = TRUE) }
      if (hele.info > 0) {hele = read.alignment(paste0(dir2,lfiles[lfiles %like% ".hele.popalign.fasta"]), "fasta", forceToLower = TRUE) }
      if (hpar.info > 0) {hpar = read.alignment(paste0(dir2,lfiles[lfiles %like% ".hpar.popalign.fasta"]), "fasta", forceToLower = TRUE) }
      if (hmel.info > 0) {hmel = read.alignment(paste0(dir2,lfiles[lfiles %like% ".hmelOut.popalign.fasta"]), "fasta", forceToLower = TRUE) }
      # calculate pairwise dnds within pops and get median
      hall.md.dnds = NaN; hall.md.dn = NaN; hall.md.ds = NaN; hallint = 0;
      trio.md.dnds = NaN; trio.md.dn = NaN; trio.md.ds = NaN; trioint = 0;
      hhec.md.dnds = NaN; hhec.md.dn = NaN; hhec.md.ds = NaN; hhecint = 0;
      hele.md.dnds = NaN; hele.md.dn = NaN; hele.md.ds = NaN; heleint = 0;
      hpar.md.dnds = NaN; hpar.md.dn = NaN; hpar.md.ds = NaN; hparint = 0;
      hmel.md.dnds = NaN; hmel.md.dn = NaN; hmel.md.ds = NaN; hmelint = 0;
      # guarantee that i) all sequences are multiples of 3 [hecint] and all have the same length [hecunq]
      if (hall.info > 0) {hallint = max((getLength(hall$seq)/3)/(floor(getLength(hall$seq)/3))); hallunq = length(unique(getLength(hall$seq))) }
      if (trio.info > 0) {trioint = max((getLength(trio$seq)/3)/(floor(getLength(trio$seq)/3))); triounq = length(unique(getLength(trio$seq))) }
      if (hhec.info > 0) {hhecint = max((getLength(hhec$seq)/3)/(floor(getLength(hhec$seq)/3))); hhecunq = length(unique(getLength(hhec$seq))) }
      if (hele.info > 0) {heleint = max((getLength(hele$seq)/3)/(floor(getLength(hele$seq)/3))); heleunq = length(unique(getLength(hele$seq))) }
      if (hpar.info > 0) {hparint = max((getLength(hpar$seq)/3)/(floor(getLength(hpar$seq)/3))); hparunq = length(unique(getLength(hpar$seq))) }
      if (hmel.info > 0) {hmelint = max((getLength(hmel$seq)/3)/(floor(getLength(hmel$seq)/3))); hmelunq = length(unique(getLength(hmel$seq))) }
      #
      dfComb=data.frame()
      if (hall.info > 0) {if (hall$nb > 1 & hallint == 1 & hallunq == 1) {hall.dnds = kaks(hall, rmgap = T); df2 = data.frame(spp=rep("hall", length(as.vector(hall.dnds$ka))), cn=rep(hall$nb, length(as.vector(hall.dnds$ka))), dn=as.vector(hall.dnds$ka), ds=as.vector(hall.dnds$ks), vda=as.vector(hall.dnds$vka), vds=as.vector(hall.dnds$vks), dist=as.vector(dist.alignment(hall))); dfComb=rbind(dfComb, df2)} }
      if (trio.info > 0) {if (trio$nb > 1 & trioint == 1 & triounq == 1) {trio.dnds = kaks(trio, rmgap = T); df2 = data.frame(spp=rep("trio", length(as.vector(trio.dnds$ka))), cn=rep(trio$nb, length(as.vector(trio.dnds$ka))), dn=as.vector(trio.dnds$ka), ds=as.vector(trio.dnds$ks), vda=as.vector(trio.dnds$vka), vds=as.vector(trio.dnds$vks), dist=as.vector(dist.alignment(trio))); dfComb=rbind(dfComb, df2)} }
      if (hhec.info > 0) {if (hhec$nb > 1 & hhecint == 1 & hhecunq == 1) {hhec.dnds = kaks(hhec, rmgap = T); df2 = data.frame(spp=rep("hhec", length(as.vector(hhec.dnds$ka))), cn=rep(hhec$nb, length(as.vector(hhec.dnds$ka))), dn=as.vector(hhec.dnds$ka), ds=as.vector(hhec.dnds$ks), vda=as.vector(hhec.dnds$vka), vds=as.vector(hhec.dnds$vks), dist=as.vector(dist.alignment(hhec))); dfComb=rbind(dfComb, df2)} }
      if (hele.info > 0) {if (hele$nb > 1 & heleint == 1 & heleunq == 1) {hele.dnds = kaks(hele, rmgap = T); df2 = data.frame(spp=rep("hele", length(as.vector(hele.dnds$ka))), cn=rep(hele$nb, length(as.vector(hele.dnds$ka))), dn=as.vector(hele.dnds$ka), ds=as.vector(hele.dnds$ks), vda=as.vector(hele.dnds$vka), vds=as.vector(hele.dnds$vks), dist=as.vector(dist.alignment(hele))); dfComb=rbind(dfComb, df2)} }
      if (hpar.info > 0) {if (hpar$nb > 1 & hparint == 1 & hparunq == 1) {hpar.dnds = kaks(hpar, rmgap = T); df2 = data.frame(spp=rep("hpar", length(as.vector(hpar.dnds$ka))), cn=rep(hpar$nb, length(as.vector(hpar.dnds$ka))), dn=as.vector(hpar.dnds$ka), ds=as.vector(hpar.dnds$ks), vda=as.vector(hpar.dnds$vka), vds=as.vector(hpar.dnds$vks), dist=as.vector(dist.alignment(hpar))); dfComb=rbind(dfComb, df2)} }
      if (hmel.info > 0) {if (hmel$nb > 1 & hmelint == 1 & hmelunq == 1) {hmel.dnds = kaks(hmel, rmgap = T); df2 = data.frame(spp=rep("hmel", length(as.vector(hmel.dnds$ka))), cn=rep(hmel$nb, length(as.vector(hmel.dnds$ka))), dn=as.vector(hmel.dnds$ka), ds=as.vector(hmel.dnds$ks), vda=as.vector(hmel.dnds$vka), vds=as.vector(hmel.dnds$vks), dist=as.vector(dist.alignment(hmel))); dfComb=rbind(dfComb, df2)} }
      #
      df1 = data.frame(gene = rep(gn, nrow(dfComb)), exon = rep(ex, nrow(dfComb)))
      dfComb = cbind(df1, dfComb)
      dfWithinPop = rbind(dfWithinPop, dfComb)
    }
  }
}
# add extra variables 
dfWithinPop$dnds = dfWithinPop$dn/dfWithinPop$ds
dfWithinPop$dnds = ifelse(dfWithinPop$dnds < 0, 0, dfWithinPop$dnds)
dfWithinPop$dnds[!is.finite(dfWithinPop$dnds)] = NA
dfWithinPop$dnds[is.na(dfWithinPop$dnds)] = -0.1
dfWithinPop$Z=(dfWithinPop$dn-dfWithinPop$ds)/sqrt(dfWithinPop$vds+dfWithinPop$vda)
dfWithinPop$geneExon = paste0(dfWithinPop$gene,".",dfWithinPop$exon)
dfWithinPop$geneExon = factor(dfWithinPop$geneExon, levels = unique(dfWithinPop$geneExon))
dfWithinPop$gene = factor(dfWithinPop$gene, levels = unique(dfWithinPop$gene))
dfWithinPop$spp = factor(dfWithinPop$spp, levels=c("hmel","hhec","hele","hpar","trio","hall"))
# remove hall and trio
dfWithinPopSub = subset(dfWithinPop, spp %in% c("hmel","hhec","hele","hpar"))

subset(dfWithinPopSub, Z > 1.96)
subset(dfWithinPopSub, Z < -1.96 & spp == "hpar")
subset(dfWithinPopSub, Z < -1.96 & spp == "hele")
subset(dfWithinPopSub, Z < -1.96 & spp == "hhec")

hist(subset(dfWithinPopSub, spp == "hhec" & dnds >= 0 & dnds <= 2)$dnds, breaks = seq(0,2,0.05))
hist(subset(dfWithinPopSub, spp == "hele" & dnds >= 0 & dnds <= 2)$dnds, breaks = seq(0,2,0.05))
hist(subset(dfWithinPopSub, spp == "hpar" & dnds >= 0 & dnds <= 2)$dnds, breaks = seq(0,2,0.05))

hhecsub = subset(dfWithinPopSub, spp == "hhec" & dnds >= 0)
helesub = subset(dfWithinPopSub, spp == "hele" & dnds >= 0)
hparsub = subset(dfWithinPopSub, spp == "hpar" & dnds >= 0)
hhecsubdnds = aggregate(hhecsub$dnds, list(hhecsub$geneExon), median)
helesubdnds = aggregate(helesub$dnds, list(helesub$geneExon), median)
hparsubdnds = aggregate(hparsub$dnds, list(hparsub$geneExon), median)

plot(hhecsubdnds$Group.1, hhecsubdnds$x)
plot(helesubdnds$Group.1, helesubdnds$x)
plot(hparsubdnds$Group.1, hparsubdnds$x)


#ggplot(dfWithinPop, aes(x=cn, y=dist, colour=spp)) + geom_point(alpha=0.10) + geom_smooth(method = "lm", col="black", alpha=0.1) + facet_wrap(~spp)

# plot
# dn
ggplot(dfWithinPop, aes(x=geneExon, y=dn, colour=gene, fill=gene)) + 
  geom_boxplot(alpha=0.1) + geom_jitter(alpha=0.5) +
  facet_wrap(~spp, nrow=5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
# ds
ggplot(dfWithinPop, aes(x=geneExon, y=ds, colour=gene, fill=gene)) + 
  geom_boxplot(alpha=0.1) + geom_jitter(alpha=0.5) +
  facet_wrap(~spp, nrow=5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
# dn/ds
ggplot(dfWithinPop, aes(x=geneExon, y=dnds, colour=gene, fill=gene)) + 
  # geom_boxplot(alpha=0.1) + 
  geom_jitter(alpha=0.5) +
  facet_wrap(~spp, nrow=3) + ylim(-0.1,0.5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
# distance
ggplot(dfWithinPop, aes(x=geneExon, y=(dist^2), colour=gene, fill=gene)) + 
  geom_boxplot(alpha=0.1) + geom_jitter(alpha=0.1) +
  facet_wrap(~spp, nrow=5) + ylab("Absolute Pairwise Distance") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

##
ggplot(dfWithinPop, aes(x=spp, y=dnds, colour=gene, fill=gene)) + 
  geom_boxplot(alpha=0.1) + geom_jitter(alpha=0.1) +
  facet_wrap(~geneExon, nrow=5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

ggplot(dfWithinPop, aes(x=spp, y=(dist^2), colour=gene, fill=gene)) + 
  geom_boxplot(alpha=0.1) + geom_jitter(alpha=0.1) +
  facet_wrap(~geneExon) + ylab("Absolute Pairwise Distance") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

# aaa = subset(dfWithinPop, gene %in% c("HMEL008474g1.t1","HMEL008474g1.t2"))
aaa = subset(dfWithinPop, gene %in% c("HMEL037926g1.t1"))
ggplot(aaa, aes(x=geneExon, y=dnds, colour=gene, fill=gene)) + 
  geom_boxplot(alpha=0.1) +
  geom_jitter(alpha=0.5) +
  facet_wrap(~spp, nrow=5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

# long
dfWithinPopLong = melt(dfWithinPop, id.vars=c("gene","exon","spp","geneExon"))
dfWithinPopLong$variable = factor(dfWithinPopLong$variable, levels = c("dn","ds","dnds","dist"))
ggplot(dfWithinPopLong, aes(x=geneExon, y=value, fill=gene, colour=gene)) + 
  geom_boxplot(alpha=0.1) +
  geom_jitter(alpha=0.25) +
  facet_wrap(~spp+variable, scales="free_y") + 
  # theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_blank())
