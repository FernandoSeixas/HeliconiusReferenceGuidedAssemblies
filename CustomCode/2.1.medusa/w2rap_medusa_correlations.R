require(ggplot2)
require(gridExtra)
options(scipen = 999)

## read data
medusaStats = read.table("/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/medusa_stats.txt", header=T)

## test correlation
cor.test(medusaStats$w2rap_contigs, medusaStats$medusa_contigs, method = "spearman")
cor.test(medusaStats$w2rap_contigs, medusaStats$medusa_lenTotal, method = "spearman")
cor.test(medusaStats$w2rap_contigs, medusaStats$medusa_lenAnchor, method = "spearman")
cor.test(medusaStats$w2rap_N50, medusaStats$medusa_contigs, method = "spearman")
cor.test(medusaStats$w2rap_N50, medusaStats$medusa_lenTotal, method = "spearman")
cor.test(medusaStats$w2rap_N50, medusaStats$medusa_lenAnchor, method = "spearman")

## plot correlation
a1 = ggplot(medusaStats, aes(x=w2rap_contigs, y=medusa_contigs)) + geom_point(col=ifelse(medusaStats$assembly %in% c("hhim","hhimfat"), "orange", "black")) + geom_smooth(method = "lm") + theme_bw() + ylim(0,5100)
a2 = ggplot(medusaStats, aes(x=w2rap_contigs, y=medusa_lenTotal/1000000)) + geom_point(col=ifelse(medusaStats$assembly %in% c("hhim","hhimfat"), "orange", "black")) + geom_smooth(method = "lm") + theme_bw() + ylim(250,450)
a3 = ggplot(medusaStats, aes(x=w2rap_contigs, y=medusa_lenAnchor/1000000)) + geom_point(col=ifelse(medusaStats$assembly %in% c("hhim","hhimfat"), "orange", "black")) + geom_smooth(method = "lm") + theme_bw() + ylim(250,450)
b1 = ggplot(medusaStats, aes(x=w2rap_N50, y=medusa_contigs)) + geom_point(col=ifelse(medusaStats$assembly %in% c("hhim","hhimfat"), "orange", "black")) + geom_smooth(method = "lm") + theme_bw() + ylim(0,5100)
b2 = ggplot(medusaStats, aes(x=w2rap_N50, y=medusa_lenTotal/1000000)) + geom_point(col=ifelse(medusaStats$assembly %in% c("hhim","hhimfat"), "orange", "black")) + geom_smooth(method = "lm") + theme_bw() + ylim(250,450)
b3 = ggplot(medusaStats, aes(x=w2rap_N50, y=medusa_lenAnchor/1000000)) + geom_point(col=ifelse(medusaStats$assembly %in% c("hhim","hhimfat"), "orange", "black")) + geom_smooth(method = "lm") + theme_bw() + ylim(250,450)

p = grid.arrange(
  a1,b1,
  a2,b2,
  a3,b3
  )
p

ggsave(filename = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/w2rap_medusa_Correlations.eps", plot = p, dpi = 300)
