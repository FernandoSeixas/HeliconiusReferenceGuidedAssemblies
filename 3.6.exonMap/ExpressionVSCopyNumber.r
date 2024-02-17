## load modules
require(stringr)
require(tidyr)
require(ggplot2)
require(ggrepel)


#################### Copy Number Variation ####################
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/1.Map2CN/0.minimap2/"
lfiles = list.files(dir)
lfiles = lfiles[grepl(x = lfiles, pattern = c("par"))]
lfiles = lfiles[grepl(x = lfiles, pattern = c("hmel02|hmel04|hmel08|hmel09"))]
dfExonMapsFilter = data.frame()
for (file in lfiles) {
  inp = read.table(paste0(dir,file), fill = TRUE)
  dfExonMapsFilter = rbind(dfExonMapsFilter, inp)
}

names(dfExonMapsFilter)[1] = "gid"
names(dfExonMapsFilter)[6] = "RGAinfo"


dfExonMapsFilter$gid = str_replace(dfExonMapsFilter$gid, "\\_RA\\-", "\\.RA\\_")
dfExonMapsFilter = dfExonMapsFilter %>% separate(gid, c("transcript_id", "exon","chrom","sta","end"), sep="_")
dfExonMapsFilter = dfExonMapsFilter %>% separate(RGAinfo, c("RGAchrom", "RGAscaffold"), sep="-")
dfExonMapsFilter$REFchromCode = substr(dfExonMapsFilter$chrom, 6, 7) 
dfExonMapsFilter$RGAchromCode = substr(dfExonMapsFilter$RGAchrom, 5, 6) 
dfExonMapsFilter = subset(dfExonMapsFilter, REFchromCode == RGAchromCode)
dfExonMapsFilter$exon_id = paste0(dfExonMapsFilter$transcript_id,".",dfExonMapsFilter$exon)
dfExonMapsFilter$Species = substr(dfExonMapsFilter$RGAchrom, 1, 4)
dfExonMapsFilter$dummy = 1


#####
CNV = aggregate(dfExonMapsFilter$dummy, list(paste0(dfExonMapsFilter$Species,".",dfExonMapsFilter$exon_id)), sum)
names(CNV) = c("id", "CN")
CNV = CNV %>% separate(id, c("Species", "gene", "transcript", "exon"), "([.])") 
CNV$transcript_id = paste0(CNV$gene,".",CNV$transcript)
CNV$exon_id = paste0(CNV$gene,".",CNV$transcript,".",CNV$exon)


#####
CNV_transcript = aggregate(CNV$CN, list(CNV$transcript_id), max)
CNV_transcript = aggregate(CNV$CN, list(CNV$transcript_id), median)
# CNV_transcript = aggregate(CNV$CN, list(CNV$transcript_id), mean)
names(CNV_transcript) = c("target_id", "CN")



######################### Expression Data #########################
ExpressionDataFile = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/nateExpressionData/duplication.30Jun20.diffExpression.txt"
ExpressionData = read.table(ExpressionDataFile, header=T)

## add CN information
ExpressionDataCN = merge(ExpressionData, CNV_transcript)

## add gene names
geneNames = read.table("/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/nateExpressionData/geneNames")
names(geneNames) = c("gene", "flybaseName")
ExpressionDataCN$gene = substr(ExpressionDataCN$target_id, 1, 12)
ExpressionDataCN = merge(ExpressionDataCN, geneNames)

## correlation test
ctest = cor.test(x=ExpressionDataCN$CN, y=ExpressionDataCN$b)
val1 = as.numeric(round(ctest$estimate, 2))
val2 = round(ctest$p.value,2)

## plot
ggplot(ExpressionDataCN, aes(x=CN, y=b, label=paste0(target_id,"\n(",flybaseName,")"))) + 
  geom_smooth(method = "lm") +
  geom_point(pch=ifelse(ExpressionDataCN$qval < 0.05, 19, 1), size=3) +
  geom_text_repel(alpha=ifelse(ExpressionDataCN$qval < 0.05, 1, 0), size=4) +
  geom_text(aes(x=1, y=5), label=paste0("Spearman's Rank Correlation Test:\nRho=",val1," (P-value=",val2,")"), 
    size=4, hjust = 0) +
  scale_x_continuous(breaks = seq(0,20,1)) +
  theme(panel.grid.minor.x = element_blank()) +
  ylab("b = log2(counts + 0.5)")
ggsave(paste0(dir,"ExpressionVsCopyNumber.eps"), dpi = 150)


