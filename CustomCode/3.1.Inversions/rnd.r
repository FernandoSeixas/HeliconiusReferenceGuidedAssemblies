require(ggplot2)
require(gridExtra)


## read data
ref = "hmelv25"
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/fdM/5.pdistance/"
fil = paste0("wgenome.",ref,".4species.dist.csv")
inp = read.csv(paste0(dir,fil))
# add dummy pos 
inp$dummypos = seq(1,nrow(inp),1)
# filter data
plot(density(inp$sites))
inpFilter = subset(inp, sites >= 10000)
nrow(inpFilter)/nrow(inp)


# get chromosome names
if (ref == "hmelv25") { inpFilter$chrom = substr(inpFilter$scaffold, 1, 7) }
if (ref == "heradem") { inpFilter$chrom = substr(inpFilter$scaffold, 1, 8) }
odds = unique(inpFilter$chrom)[seq(1,21,2)]
inpFilter$chromCol = ifelse(inpFilter$chrom %in% odds, "black", "darkgrey")

# determine chromosome middle positions
chromLabels = data.frame(chrom = unique(inpFilter$chrom), midpos = 0)
for (i in 1:length(unique(inpFilter$chrom))) {
  ch = unique(inpFilter$chrom)[i]
  s1 = subset(inpFilter, chrom == as.character(ch))
  chromLabels$chrom[i] = ch
  chromLabels$midpos[i] = median(s1$dummypos)
}

# rnd
inpFilter$RND_bur_era_mel = inpFilter$dxy_era_bur/inpFilter$dxy_era_mel
inpFilter$RND_dor_era_mel = inpFilter$dxy_era_dor/inpFilter$dxy_era_mel

# calculate RND z.scores
inpFilter$RND_bur_era_mel.zscore = (inpFilter$RND_bur_era_mel-mean(inpFilter$RND_bur_era_mel)) / sd(inpFilter$RND_bur_era_mel)
inpFilter$RND_dor_era_mel.zscore = (inpFilter$RND_dor_era_mel-mean(inpFilter$RND_dor_era_mel)) / sd(inpFilter$RND_dor_era_mel)

# zscore limits [for plot]
zscore.low.bur = mean(inpFilter$RND_bur_era_mel) - 3*sd(inpFilter$RND_bur_era_mel)
zscore.hig.bur = mean(inpFilter$RND_bur_era_mel) + 3*sd(inpFilter$RND_bur_era_mel)
zscore.low.dor = mean(inpFilter$RND_dor_era_mel) - 3*sd(inpFilter$RND_dor_era_mel)
zscore.hig.dor = mean(inpFilter$RND_dor_era_mel) + 3*sd(inpFilter$RND_dor_era_mel)

# breakpoints
if (ref == "heradem") {
  chr2 = subset(inp, scaffold == "Herato0206" & start < 1059137 & end > 880283)
  chr6 = subset(inp, scaffold == "Herato0606" & start < 12756548 & end > 12189406)
  chr13 = subset(inp, scaffold == "Herato1301" & start < 22668029 & end > 22307500)
  chr21 = subset(inp, scaffold == "Herato2101" & start < 10815508 & end > 10443345)
}
if (ref == "hmelv25") {
  chr2 = subset(inp, scaffold == "Hmel202001o" & start < 2230606 & end > 2147755)
  chr6 = subset(inp, scaffold == "Hmel206001o" & start < 12591271 & end > 12129405)
  chr13 = subset(inp, scaffold == "Hmel213001o" & start < 17154514 & end > 16917575)
  chr21 = subset(inp, scaffold == "Hmel221001o" & start < 8342578 & end > 8021528)
}
# plot RND
p1 <- ggplot(inpFilter) + 
  geom_hline(yintercept = c(zscore.low.bur,zscore.hig.bur), col="red", alpha=0.5) +
  geom_point(aes(x=dummypos, y=RND_bur_era_mel), alpha=ifelse(abs(inpFilter$RND_bur_era_mel.zscore) >= 3, 1.0, 0.4), col=inpFilter$chromCol) +
  scale_x_continuous(breaks = chromLabels$midpos, labels = seq(1,21,1)) +
  geom_vline(xintercept = c(min(chr2$dummypos), max(chr2$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr6$dummypos), max(chr6$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr13$dummypos), max(chr13$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr21$dummypos), max(chr21$dummypos)), alpha=0.2, lty=2) +
  xlab("Chromosome") + ylab ("RND") + ylim(0.5,2.2) + ggtitle("C") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 9, angle=90, hjust=1.0, vjust=0.5),
        axis.text.y = element_text(size = 9, angle=0, hjust=1.0, vjust=0.5),
        axis.title = element_text(size = 14)
  )
p2 <- ggplot(inpFilter) + 
  geom_hline(yintercept = c(zscore.low.dor,zscore.hig.dor), col="red", alpha=0.5) +
  geom_point(aes(x=dummypos, y=RND_dor_era_mel), alpha=ifelse(abs(inpFilter$RND_dor_era_mel.zscore) >= 3, 1.0, 0.4), col=inpFilter$chromCol) +
  scale_x_continuous(breaks = chromLabels$midpos, labels = seq(1,21,1)) +
  geom_vline(xintercept = c(min(chr2$dummypos), max(chr2$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr6$dummypos), max(chr6$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr13$dummypos), max(chr13$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr21$dummypos), max(chr21$dummypos)), alpha=0.2, lty=2) +
  xlab("Chromosome") + ylab ("RND") + ylim(0.5,2.2) + ggtitle("D") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 9, angle=90, hjust=1.0, vjust=0.5),
        axis.text.y = element_text(size = 9, angle=0, hjust=1.0, vjust=0.5),
        axis.title = element_text(size = 14)
  )
grid.arrange(p1, p2, p3, p4, nrow=2)


subset(inpFilter, RND_bur_era_mel.zscore < -3)[,1:3]
subset(inpFilter, RND_dor_era_mel.zscore < -3)[,1:3]
subset(inpFilter, RND_bur_era_mel.zscore < -3 & RND_dor_era_mel.zscore < -3)


# mmm = subset(inpFilter, scaffold == "Herato1301")
# mmm = subset(inpFilter, scaffold == "Hmel213001o")
# ggplot(mmm) + 
#   geom_line(aes(x=dummypos, y=RND_bur_era_mel+2)) +
#   geom_line(aes(x=dummypos, y=RND_dor_era_mel+0)) +
#   geom_vline(xintercept = c(min(chr13$dummypos), max(chr13$dummypos)), alpha=0.2, lty=2) +
#   ylim(0,4)

