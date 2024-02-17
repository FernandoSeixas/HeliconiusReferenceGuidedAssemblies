require(seqinr)
require(stringr)
require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)

## All Pops Dn/Ds  ==================================================
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/"
dir1 = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/1.Map2CN/5.frmseq/"
chrom = "hmel20"
exonsTable = read.table(paste0(dir,chrom,".target.exonsFrame.txt"), header=T)
#
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
        all.kaks = kaks(all, rmgap = T);
        all.dn=all.kaks$ka[1:(all$nb-1)]
        all.ds=all.kaks$ks[1:(all$nb-1)]
        all.vdn=all.kaks$vka[1:(all$nb-1)]
        all.vds=all.kaks$vks[1:(all$nb-1)]
        all.nm=all$nam[2:all$nb]
        #
        # all.nm = substr(all.nm, 1, 4)
        #
        df1 = data.frame(gene = rep(gn, (all$nb-1)), exon = rep(ex, (all$nb-1)), name=all.nm, dn = all.dn, ds = all.ds, vdn=all.vdn, vds=all.vds)
        dfAllPop = rbind(dfAllPop, df1)
      }
      if (all$nb == 1 | allint != 1) {
        df2 = data.frame(gene = gn, exon = ex, name="hmel", dn=NA, ds=NA, vdn=NA, vds=NA); dfAllPop = rbind(dfAllPop, df2)
        df2 = data.frame(gene = gn, exon = ex, name="hbes", dn=NA, ds=NA, vdn=NA, vds=NA); dfAllPop = rbind(dfAllPop, df2)
        df2 = data.frame(gene = gn, exon = ex, name="hnum", dn=NA, ds=NA, vdn=NA, vds=NA); dfAllPop = rbind(dfAllPop, df2)
        df2 = data.frame(gene = gn, exon = ex, name="hhec", dn=NA, ds=NA, vdn=NA, vds=NA); dfAllPop = rbind(dfAllPop, df2)
        df2 = data.frame(gene = gn, exon = ex, name="hele", dn=NA, ds=NA, vdn=NA, vds=NA); dfAllPop = rbind(dfAllPop, df2)
        df2 = data.frame(gene = gn, exon = ex, name="hpar", dn=NA, ds=NA, vdn=NA, vds=NA); dfAllPop = rbind(dfAllPop, df2)
        
      }
    }
  }
}

# create geneExon names
dfAllPop$gene = factor(dfAllPop$gene, levels = unique(dfAllPop$gene))
dfAllPop$geneExon = paste0(dfAllPop$gene,dfAllPop$exon)
dfAllPop$geneExon = factor(dfAllPop$geneExon, levels = unique(dfAllPop$geneExon))
#
tail(dfAllPop)
dfAllPop$name = str_remove(dfAllPop$name, "_R_")
# calculate dnds, transform negative dn/ds values to 0 and infinite values to NA
dfAllPop$dnds = dfAllPop$dn/dfAllPop$ds
dfAllPop$dnds = ifelse(dfAllPop$dnds < 0, NA, dfAllPop$dnds)
dfAllPop$dnds[!is.finite(dfAllPop$dnds)] = NA
dfAllPop$dnds[is.na(dfAllPop$dnds)] = -0.1
dfAllPop$Z = (dfAllPop$dn-dfAllPop$ds)/sqrt(dfAllPop$vdn+dfAllPop$vds)
# improve species names
dfAllPop <- dfAllPop %>% separate(name, c("spp","Ngene","transcript","Nexon","copy"),"[.]")
dfAllPop$spp = substr(dfAllPop$spp, 1, 4)
dfAllPop$newname = paste0(dfAllPop$spp,".",dfAllPop$copy)
dfAllPop$newname = str_replace(dfAllPop$spp, ".NA", ".copy1")
# dfAllPop$spp = str_remove(dfAllPop$spp, cnb)
dfAllPop$spp = factor(dfAllPop$spp, levels = c("hmel","htim","hcyd","hbes","hnum","hhec","hele","hpar","hera"))
dfAllPop = dfAllPop[dfAllPop$spp != "hera",]
# test selection hypothesis
dfAllPop$NonNeutral = ifelse(abs(dfAllPop$Z) > 1.96, "yes", "no")
dfAllPop$PurSelect = ifelse(dfAllPop$Z < -1.96, "yes", "no")
dfAllPop$PosSelect = ifelse(dfAllPop$Z >  1.96, "yes", "no")


######################### PLOT DN/DS #########################
# dnds per exon 
mx = max(dfAllPop$dnds)
ggplot(dfAllPop, aes(x=geneExon, y=dnds, fill=gene, colour=gene)) +
  geom_boxplot(alpha=0.1) + geom_jitter(alpha=0.5) +
  facet_wrap(~spp, nrow=8) + ylab("Dn/Ds\n(to hmelv25 reference sequences)") + 
  # ylim(-0.1,mx) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5), legend.position = "none") + 
  ggtitle(chrom)

# dnds per gene
dfAllPop$dnds = ifelse(dfAllPop$dnds < 0, 0, dfAllPop$dnds)
ggplot(dfAllPop, aes(x=gene, y=dnds, fill=gene, colour=gene)) +
  geom_boxplot(alpha=0.1) + 
  # geom_jitter(alpha=0.3) +
  facet_wrap(~spp, nrow=4) + 
  xlab("Gene Transcript") +
  ylab("Dn/Ds\n(to hmelv25 reference sequences)") + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5), legend.position = "none") + 
  ggtitle(chrom)




medianSppExon = aggregate(dfAllPop$dnds, list(paste0(dfAllPop$newname,".",dfAllPop$geneExon)), median)
names(medianSppExon) = c("sppExon","dnds")
medianSppExon <- separate(medianSppExon, "sppExon", c("spp","Gene","transExon"), sep = "\\.")
medianSppExon$GeneExon = paste0(medianSppExon$Gene, medianSppExon$transExon)
medianSppExon$dnds = ifelse(medianSppExon$dnds < 0, 0, medianSppExon$dnds)
head(medianSppExon)

ggplot(medianSppExon) + 
  geom_boxplot(aes(x=Gene, y=dnds, col=Gene)) + 
  facet_wrap(~spp, nrow=4) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5), legend.position = "none")
  

# dn/ds distribution across all exons
# density plot
a = subset(dfAllPop, spp %in% c("hhec","hele","hpar","hbes","hnum","hcyd","htim"))
ggplot(a) + 
  geom_histogram(aes(x=dnds, fill=spp), alpha=0.5, breaks=seq(0,mx,0.05), col="black") + 
  facet_wrap(~spp, nrow=5, scales = "free_y") + ggtitle(chrom)
# violin plot
ggplot(a) + geom_violin(aes(x=spp, y=dnds, fill=spp), alpha=0.5, col="black") + ggtitle(chrom)
ggplot(a) + geom_boxplot(aes(x=spp, y=dnds, fill=spp), alpha=0.5, col="black") + ggtitle(chrom)
# mx = max(dfAllPop$dnds, na.rm=T)
# hist(subset(dfAllPop, spp == "hhec" & dnds >= 0)$dnds, breaks=seq(0,mx+0.1,0.05), col="orange")
# hist(subset(dfAllPop, spp == "hele" & dnds >= 0)$dnds, breaks=seq(0,mx+0.1,0.05), col="blue")
# hist(subset(dfAllPop, spp == "hpar" & dnds >= 0)$dnds, breaks=seq(0,mx+0.1,0.05), col="red")
# hist(subset(dfAllPop, spp == "hbes" & dnds >= 0)$dnds, breaks=seq(0,mx+0.1,0.05), col="black")
# hist(subset(dfAllPop, spp == "hnum" & dnds >= 0)$dnds, breaks=seq(0,mx+0.1,0.05), col="grey")

plot(density(subset(dfAllPop, spp == "hhec" & dnds >= 0)$dnds, na.rm=T), col="orange", xlab="Dn/Ds", ylim=c(0,2), main=chrom)
lines(density(subset(dfAllPop, spp == "hele" & dnds >= 0)$dnds, na.rm=T), col="blue")
lines(density(subset(dfAllPop, spp == "hpar" & dnds >= 0)$dnds, na.rm=T), col="red")
lines(density(subset(dfAllPop, spp == "hbes" & dnds >= 0)$dnds, na.rm=T), col="black")
lines(density(subset(dfAllPop, spp == "hnum" & dnds >= 0)$dnds, na.rm=T), col="grey")



######################### DN/DS vs CNV #########################
head(dfAllPop)
cnvTable = data.frame(table(paste0(dfAllPop$geneExon,".",dfAllPop$spp)))
names(cnvTable) = c("geneSpp", "cnv")
dfAllPop$geneSpp = paste0(dfAllPop$geneExon,".",dfAllPop$spp)
dfAllPop$cnv = 0
# 
for (rw in 1:nrow(dfAllPop)) {
  gspp=dfAllPop$geneSpp[rw]
  dfAllPop$cnv[rw] = cnvTable[cnvTable$geneSpp == gspp,]$cnv
}
# plot
ggplot(dfAllPop, aes(x=cnv, y=dnds, colour=spp, fill=spp, group=cnv)) +
  geom_boxplot(alpha=0.05) + geom_jitter(alpha=0.5) + 
  scale_x_continuous(breaks=seq(0,25,1)) + facet_wrap(~spp, nrow=5) +
  geom_hline(yintercept = 1, lty=2)
  
