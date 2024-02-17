## load modules
require(tidyr)
require(reshape2)
require(stringr)
require(ggplot2)


## read data
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.3.GenomeExpansion/exonMapping/1.Map2CN/"
fil = "exonsSummary.txt"
exonSummary = read.table(paste0(dir,fil))
names(exonSummary) = 
  c("exonComplex",
    "hmel.cnv","hhec.cnv","hele.cnv","hpar.cnv",
    "hmel.stp","hhec.stp","hele.stp","hpar.stp",
    "hmel.frm","hhec.frm","hele.frm","hpar.frm",
    "hmel.sub","hhec.sub","hele.sub","hpar.sub"
  )
# split columns when necessary
exonSummary$exonComplex = str_replace(exonSummary$exonComplex, "-", ".")
exonSummary <- exonSummary %>% separate(exonComplex, c("chrom","gene","transcript","exon"), sep="\\.")
exonSummary$geneExon = paste0(exonSummary$gene,".",exonSummary$transcript,".",exonSummary$exon)


######################### Analyze the data #########################
# split by type [cnv, stp, frm, sub]
cnv <- melt(exonSummary[,c(1:4, 21, 5:8 )], names(exonSummary[c(1:4,21)])); names(cnv)[6] = "spp"; names(cnv)[7] = "cnv"; cnv$spp = str_remove(cnv$spp, ".cnv")
stp <- melt(exonSummary[,c(1:4, 21, 9:12)], names(exonSummary[c(1:4,21)])); names(stp)[6] = "spp"; names(stp)[7] = "stp"; stp$spp = str_remove(stp$spp, ".stp")
frm <- melt(exonSummary[,c(1:4, 21,13:16)], names(exonSummary[c(1:4,21)])); names(frm)[6] = "spp"; names(frm)[7] = "frm"; frm$spp = str_remove(frm$spp, ".frm")
sub <- melt(exonSummary[,c(1:4, 21,17:20)], names(exonSummary[c(1:4,21)])); names(sub)[6] = "spp"; names(sub)[7] = "sub"; sub$spp = str_remove(sub$spp, ".sub")
# combine dataframe
cmb = merge(cnv, stp, by = names(cnv)[1:6])
cmb = merge(cmb, frm, by = names(cnv)[1:6])
cmb = merge(cmb, sub, by = names(cnv)[1:6])
# add transcript_id
cmb$transcript_id = paste0(cmb$gene,".",cmb$transcript)
# determine species
cmb$spp = factor(cmb$spp, levels = c("hmel","hhec","hele","hpar"))



## Proportion of exons with stop codons
# Proportion of exons
nrow(subset(cmb, spp == "hmel" & stp > 0))/nrow(subset(cmb, spp == "hmel"))*100
nrow(subset(cmb, spp == "hhec" & stp > 0))/nrow(subset(cmb, spp == "hhec"))*100
nrow(subset(cmb, spp == "hele" & stp > 0))/nrow(subset(cmb, spp == "hele"))*100
nrow(subset(cmb, spp == "hpar" & stp > 0))/nrow(subset(cmb, spp == "hpar"))*100
# Proportion of exon copies 
round(sum(subset(cmb, spp == "hmel")$stp)/sum(subset(cmb, spp == "hmel")$cnv)*100)
round(sum(subset(cmb, spp == "hhec")$stp)/sum(subset(cmb, spp == "hhec")$cnv)*100)
round(sum(subset(cmb, spp == "hele")$stp)/sum(subset(cmb, spp == "hele")$cnv)*100)
round(sum(subset(cmb, spp == "hpar")$stp)/sum(subset(cmb, spp == "hpar")$cnv)*100)
# plot
# ggplot(cmb) + 
#   # geom_point(aes(x=geneExon, y=stp, colour=spp), alpha=0.2) +
#   geom_boxplot(aes(x=transcript_id, y=stp/cnv, fill=spp, colour=spp), alpha=0.2) +
#   facet_wrap(~spp, nrow=4) +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))


## per chromosome /////
chr = "hmel09"
# get exon names to then order in plot
geneExonNames = read.table(paste0(dir,"../",chr,".target.geneExons.txt"))[,1]
gene_levels = unique(sapply(str_split(geneExonNames,"\\."), `[`, 1))
gene_levels = unique(sapply(str_split(gene_levels,"-"), `[`, 1))
# subset data to specific chromosome and reorder genes
cmbSub = subset(cmb, chrom == chr)
cmbSubOrder = data.frame()
for (gid in gene_levels) {
  s1 = subset(cmbSub, gene == gid)
  cmbSubOrder = rbind(cmbSubOrder, s1)
}
cmbSubOrder$gene = factor(cmbSubOrder$gene, levels = unique(cmbSubOrder$gene))
cmbSubOrder$geneExon = factor(cmbSubOrder$geneExon, levels = unique(cmbSubOrder$geneExon))
cmbSubOrder$transcript_id = str_remove(cmbSubOrder$transcript_id, "\\.RA")
cmbSubOrder$transcript_id = factor(cmbSubOrder$transcript_id, levels = unique(cmbSubOrder$transcript_id))
# plot
ggplot(cmbSubOrder) + 
  # geom_boxplot(aes(x=transcript_id, y=stp/cnv, fill=transcript_id, colour=transcript_id), alpha=0.2) +
  # ylim(0,1) +
  geom_boxplot(aes(x=geneExon, y=cnv), alpha=0.2) +
  geom_col(aes(x=geneExon, y=cnv-stp, fill=transcript_id, colour=transcript_id), alpha=0.2) +
  geom_hline(yintercept = 1, lty=2) +
  facet_wrap(~spp, nrow=4) + 
  theme(axis.text.x = element_blank()) +
  # theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  xlab("Exon") + ylab("Copy Number (CN)") + 
  labs(fill="Gene Transcript", col="Gene Transcript")




# ggplot(cmbSubOrder) + 
#   # geom_boxplot(aes(x=transcript_id, y=stp/cnv, fill=transcript_id, colour=transcript_id), alpha=0.2) +
#   # ylim(0,1) +
#   geom_boxplot(aes(x=transcript_id, y=cnv-stp, fill=transcript_id, colour=transcript_id), alpha=0.1) +
#   geom_boxplot(aes(x=transcript_id, y=cnv), alpha=0.1) +
#   # geom_hline(yintercept = 1, lty=2) +
#   facet_wrap(~spp, nrow=4) + 
#   theme(axis.text.x = element_blank()) +
#   # theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
#   xlab("Gene Transcript") + ylab("Copy Number (CN)") + 
#   labs(fill="Gene Transcript", col="Gene Transcript")
# 
# ggplot(cmbSubOrder) + 
#   geom_point(aes(x=geneExon, y=cnv), color="black", alpha=0.6) +
#   geom_point(aes(x=geneExon, y=stp, colour=transcript_id), alpha=0.6) +
#   # geom_point(aes(x=geneExon, y=cnv-stp, colour=transcript_id), alpha=0.8) +
#   geom_hline(yintercept = 1, lty=2) +
#   facet_wrap(~spp, nrow=4, scales = "free_y") +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
#   xlab("Exon") + ylab("Proportion of Exon Copies with Stop Codons")