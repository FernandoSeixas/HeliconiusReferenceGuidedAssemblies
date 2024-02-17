## load modules
require(ggplot2)
require(gridExtra)

## read data  ==================================================
ref = "hmelv25"
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/fdM/4.dstats/"
fil1 = paste0("wgenome.",ref,".outgroup.helbur.csv")
fil2 = paste0("wgenome.",ref,".outgroup.lapdor.csv")
inp1 = read.csv(paste0(dir,fil1))
inp2 = read.csv(paste0(dir,fil2))
inpFilter1 = subset(inp1, sitesUsed >= 200)
inpFilter2 = subset(inp2, sitesUsed >= 200)

# combine datasets
inpFilterCombined = merge(inpFilter1, inpFilter2, by = c("scaffold", "start", "end"))

# reorder windows by position
inpFilterCombinedReorder = data.frame()
for (sc in unique(inpFilterCombined$scaffold)) {
  s1 = subset(inpFilterCombined, scaffold == sc)
  s1 = s1[order(s1$start),]
  inpFilterCombinedReorder = rbind(inpFilterCombinedReorder, s1)
}
inpFilterCombined = inpFilterCombinedReorder
inpFilterCombined$dummypos = seq(1, nrow(inpFilterCombined), 1)

# calculate fdm and fd Zscores
inpFilterCombined$fdM.x.zscore = (inpFilterCombined$fdM.x-mean(inpFilterCombined$fdM.x)) / sd(inpFilterCombined$fdM.x)
inpFilterCombined$fdM.y.zscore = (inpFilterCombined$fdM.y-mean(inpFilterCombined$fdM.y)) / sd(inpFilterCombined$fdM.y)
inpFilterCombined$fd.x.zscore = (inpFilterCombined$fd.x-mean(inpFilterCombined$fd.x)) / sd(inpFilterCombined$fd.x)
inpFilterCombined$fd.y.zscore = (inpFilterCombined$fd.y-mean(inpFilterCombined$fd.y)) / sd(inpFilterCombined$fd.y)

# get chromosome names
if (ref == "hmelv25") { inpFilterCombined$chrom = substr(inpFilterCombined$scaffold, 1, 7) }
if (ref == "heradem") { inpFilterCombined$chrom = substr(inpFilterCombined$scaffold, 1, 8) }
odds = unique(inpFilterCombined$chrom)[seq(1,21,2)]
inpFilterCombined$chromCol = ifelse(inpFilterCombined$chrom %in% odds, "black", "darkgrey")

# determine chromosome middle positions
chromLabels = data.frame(chrom = unique(inpFilterCombined$chrom), midpos = 0)
for (i in 1:length(unique(inpFilterCombined$chrom))) {
  ch = unique(inpFilterCombined$chrom)[i]
  s1 = subset(inpFilterCombined, chrom == as.character(ch))
  chromLabels$chrom[i] = ch
  chromLabels$midpos[i] = median(s1$dummypos)
}




## plot  ===============================================================
zscore.low.x = mean(inpFilterCombined$fdM.x) - 3*sd(inpFilterCombined$fdM.x)
zscore.hig.x = mean(inpFilterCombined$fdM.x) + 3*sd(inpFilterCombined$fdM.x)
zscore.low.y = mean(inpFilterCombined$fdM.y) - 3*sd(inpFilterCombined$fdM.y)
zscore.hig.y = mean(inpFilterCombined$fdM.y) + 3*sd(inpFilterCombined$fdM.y)

# # highlight zone
# hl = subset(inpFilterCombined, fdM.x.zscore < -3 & fdM.y.zscore < -3)

# breakpoints
if (ref == "heradem") {
  chr2 = subset(inpFilterCombined, scaffold == "Herato0206" & start < 1059137 & end > 880283)
  chr6 = subset(inpFilterCombined, scaffold == "Herato0606" & start < 12756548 & end > 12189406)
  chr13 = subset(inpFilterCombined, scaffold == "Herato1301" & start < 22668029 & end > 22307500)
  chr21 = subset(inpFilterCombined, scaffold == "Herato2101" & start < 10815508 & end > 10443345)
}
if (ref == "hmelv25") {
  chr2 = subset(inpFilterCombined, scaffold == "Hmel202001o" & start < 2230606 & end > 2147755)
  chr6 = subset(inpFilterCombined, scaffold == "Hmel206001o" & start < 12591271 & end > 12129405)
  chr13 = subset(inpFilterCombined, scaffold == "Hmel213001o" & start < 17154514 & end > 16917575)
  chr21 = subset(inpFilterCombined, scaffold == "Hmel221001o" & start < 8342578 & end > 8021528)
}


# plot
p3 <- ggplot(inpFilterCombined) + 
  geom_hline(yintercept = c(zscore.low.x,zscore.hig.x), col="red", alpha=0.5) +
  geom_point(aes(x=dummypos, y=fdM.x+0), alpha=ifelse(abs(inpFilterCombined$fdM.x.zscore) >= 3, 1.0, 0.4), col=inpFilterCombined$chromCol) +
  geom_vline(xintercept = c(min(chr2$dummypos), max(chr2$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr6$dummypos), max(chr6$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr13$dummypos), max(chr13$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr21$dummypos), max(chr21$dummypos)), alpha=0.2, lty=2) +
  scale_x_continuous(breaks = chromLabels$midpos, labels = seq(1,21,1)) +
  ylim(-0.5,0.5) + xlab("Chromosome") + ylab ("fdm") + ggtitle("C") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 9, angle=90, hjust=1.0, vjust=0.5),
        axis.text.y = element_text(size = 9, angle=0, hjust=1.0, vjust=0.5),
        axis.title = element_text(size = 14)
  )

p4 <- ggplot(inpFilterCombined) + 
  geom_hline(yintercept = c(zscore.low.y,zscore.hig.y), col="red", alpha=0.5) +
  geom_point(aes(x=dummypos, y=fdM.y), alpha=ifelse(abs(inpFilterCombined$fdM.y.zscore) >= 3, 1.0, 0.4), col=inpFilterCombined$chromCol) +
  geom_vline(xintercept = c(min(chr2$dummypos), max(chr2$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr6$dummypos), max(chr6$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr13$dummypos), max(chr13$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr21$dummypos), max(chr21$dummypos)), alpha=0.2, lty=2) +
  scale_x_continuous(breaks = chromLabels$midpos, labels = seq(1,21,1)) +
  ylim(-0.5,0.5) + xlab("Chromosome") + ylab ("fdm") + ggtitle("D") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 9, angle=90, hjust=1.0, vjust=0.5),
        axis.text.y = element_text(size = 9, angle=0, hjust=1.0, vjust=0.5),
        axis.title = element_text(size = 14)
  )
a = grid.arrange(p3, p1, p4, p2, ncol=2)
ggsave(file = paste0(dir,"fdm_and_rnd.svg"), plot = a, dpi = "print")


# get outliers
subset(inpFilterCombined, fdM.x.zscore > 3 & fdM.y.zscore > 3)
subset(inpFilterCombined, fdM.x.zscore < -3 & fdM.y.zscore < -3)
median(inpFilterCombined$fdM.x)
median(inpFilterCombined$fdM.y)



#################### Significance Test ####################

require("bootstrap")

i = rnorm(100)
mmm = jackknife(i, median)
plot(density(i))
plot(density(mmm$jack.values))

median(i)


bur.jk = jackknife(inpFilterCombined$fdM.x, mean)
dor.jk = jackknife(inpFilterCombined$fdM.y, mean)

plot(density(inpFilterCombined$fdM.x))

plot(density(bur.jk$jack.values))
plot(density(dor.jk$jack.values))


# ## plot specific scaffold ========================================
# # chr = "Herato0206"; bk1 = 880283; bk2 = 1059137
# # chr = "Herato0606"; bk1 = 12189406; bk2 = 12756548; bk3 = 12361782; bk4 = 12367112
# # chr = "Herato1301"; bk1 = 22307500; bk2 = 22668029
# # chr = "Herato2101"; bk1 = 10443345; bk2 = 10815508
# 
# chr="Hmel202001o"; bk1=2147755; bk2=2230606
# # chr="Hmel206001o"; bk1=12129405; bk2=12591271; bk3=12279100; bk4=12280171
# # chr="Hmel213001o"; bk1=16917575; bk2=17154514
# # chr="Hmel221001o"; bk1=8021528; bk2=8342578
# sub = subset(inpFilterCombined, scaffold == chr)
# ggplot(sub) + 
#   geom_vline(xintercept = c(bk1,bk2)/1000000, lty=2, alpha=0.5) +
#   # geom_vline(xintercept = c(bk3,bk4)/1000000, lty=2, alpha=0.5) +
#   geom_point(aes(x=(start+end)/2000000, y=fdM.x), alpha=0.4, col="blue") +
#   geom_point(aes(x=(start+end)/2000000, y=fdM.y), alpha=0.4, col="red") +
#   geom_line(aes(x=(start+end)/2000000, y=fdM.x), alpha=0.5, col="blue") +
#   geom_line(aes(x=(start+end)/2000000, y=fdM.y), alpha=0.5, col="red") +
#   scale_x_continuous(limits = c(0,max(sub$end)/1000000), breaks = seq(0,24,1), expand = c(0, 0)) +
#   ylim(-0.5,0.5) + xlab("Scaffold Position (Mb)") + ylab("fdM") +
#   theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
#   ggtitle(chr)
# 
# 
# 
# 
# aaa = subset(inpFilterCombined, fdM.x.zscore < -3 & fdM.y.zscore < -3 & scaffold == "Herato1301")
# bbb = inpFilterCombined
# ggplot() +
#   geom_boxplot(data=aaa, aes(x=1, y=sitesUsed.x)) +
#   geom_boxplot(data=bbb, aes(x=2, y=sitesUsed.x))
# 
# 
# 
