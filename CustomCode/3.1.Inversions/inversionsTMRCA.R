#################### Load Modules ###############
require(stringr)
require(reshape2)
require(ggplot2)
require(gridExtra)


#################### Get Data ####################
## read data
ref = "hmelv25"
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/fdM/5.pdistance/"
fil = paste0("wgenome.",ref,".6species.dist.csv")
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



#################### Relative Distances ####################
# # within clades
# inpFilter$RND_mel_num = inpFilter$dxy_mel_num/inpFilter$dxy_era_mel
# inpFilter$RND_era_sar = inpFilter$dxy_era_sar/inpFilter$dxy_era_mel
# inpFilter$RND_bur_dor = inpFilter$dxy_bur_dor/inpFilter$dxy_era_mel
# # between mid and era/mel
# inpFilter$RND_bur_era = inpFilter$dxy_era_bur/inpFilter$dxy_era_mel
# inpFilter$RND_bur_mel = inpFilter$dxy_mel_bur/inpFilter$dxy_era_mel
# inpFilter$RND_dor_era = inpFilter$dxy_era_dor/inpFilter$dxy_era_mel
# inpFilter$RND_dor_mel = inpFilter$dxy_mel_dor/inpFilter$dxy_era_mel

# within clades
inpFilter$RND_mel_num = inpFilter$dxy_mel_num/inpFilter$dxy_era_mel
inpFilter$RND_era_sar = inpFilter$dxy_era_sar/inpFilter$dxy_era_mel
inpFilter$RND_bur_dor = inpFilter$dxy_bur_dor/inpFilter$dxy_era_mel
# inpFilter$RND_era_mel = inpFilter$dxy_era_mel/inpFilter$dxy_era_mel
inpFilter$RND_sar_num = inpFilter$dxy_sar_num/inpFilter$dxy_era_mel
inpFilter$RND_era_num = inpFilter$dxy_era_num/inpFilter$dxy_era_mel
inpFilter$RND_sar_mel = inpFilter$dxy_sar_mel/inpFilter$dxy_era_mel
# between mid and era/mel
inpFilter$RND_bur_era = inpFilter$dxy_era_bur/inpFilter$dxy_era_mel
inpFilter$RND_bur_mel = inpFilter$dxy_mel_bur/inpFilter$dxy_era_mel
inpFilter$RND_dor_era = inpFilter$dxy_era_dor/inpFilter$dxy_era_mel
inpFilter$RND_dor_mel = inpFilter$dxy_mel_dor/inpFilter$dxy_era_mel

# breakpoints
if (ref == "heradem") {
  chr2 = subset(inpFilter, scaffold == "Herato0206" & start < 1059137 & end > 880283)
  chr6 = subset(inpFilter, scaffold == "Herato0606" & start < 12756548 & end > 12189406)
  chr13 = subset(inpFilter, scaffold == "Herato1301" & start < 22668029 & end > 22307500)
  chr21 = subset(inpFilter, scaffold == "Herato2101" & start < 10815508 & end > 10443345)
}
if (ref == "hmelv25") {
  chr2 = subset(inpFilter, scaffold == "Hmel202001o" & start < 2230606 & end > 2147755)
  chr6 = subset(inpFilter, scaffold == "Hmel206001o" & start < 12591271 & end > 12129405)
  chr13 = subset(inpFilter, scaffold == "Hmel213001o" & start < 17154514 & end > 16917575)
  chr21 = subset(inpFilter, scaffold == "Hmel221001o" & start < 8342578 & end > 8021528)
}
# if (ref == "hmelv25") {
#   chr2 = subset(inpFilter, scaffold == "Hmel202001o" &  start > 2147755  & end < 2230606)
#   chr6 = subset(inpFilter, scaffold == "Hmel206001o" &  start > 12129405 & end < 12591271)
#   chr13 = subset(inpFilter, scaffold == "Hmel213001o" & start > 16917575 & end < 17154514)
#   chr21 = subset(inpFilter, scaffold == "Hmel221001o" & start > 8021528  & end < 8342578)
# }



##### PLOT EACH INVERSION VERSUS GLOBAL TMRCAS #####
names(inpFilter)
plotInpFilter = inpFilter[,c(39,40,42,45,46)]
plotInpFilter = melt(plotInpFilter)
names(plotInpFilter) = c("SpeciesPair", "RND")
unique(plotInpFilter$SpeciesPair)
plotInpFilter$SpeciesPair = droplevels(plotInpFilter$SpeciesPair)
plotInpFilter$SpeciesPair = str_remove(plotInpFilter$SpeciesPair, "RND_")
plotInpFilter$SpeciesPair = str_replace(plotInpFilter$SpeciesPair, "_", "-")

xvals = c("bur-era","sar-num","bur-mel","era-sar","mel-num")
plotInpFilter$SpeciesPair = factor(plotInpFilter$SpeciesPair, levels = xvals)

unique(plotInpFilter$SpeciesPair)

p1 = ggplot(plotInpFilter) +
  # TMRCA global
  geom_violin(aes(x=SpeciesPair, y=RND), alpha=0.1) +
  # TMRCA OF INVERSIONS
  geom_point(aes(x="bur-era", y=median(chr2$RND_bur_era)), alpha=0.6, pch=4, col="red") +
  geom_point(aes(x="sar-num", y=median(chr2$RND_sar_num)), alpha=0.6, pch=4, col="red") +
  geom_point(aes(x="bur-mel", y=median(chr2$RND_bur_mel)), alpha=0.6, pch=4, col="red") +
  geom_point(aes(x="era-sar", y=median(chr2$RND_era_sar)), alpha=0.6, pch=4, col="red") +
  geom_point(aes(x="mel-num", y=median(chr2$RND_mel_num)), alpha=0.6, pch=4, col="red") +
  # TMRCA OF INVERSIONS - LINES
  geom_hline(aes(yintercept=median(chr2$RND_bur_mel)), alpha=0.6, lty=2, col="red") +
  geom_hline(aes(yintercept=median(chr2$RND_era_sar)), alpha=0.6, lty=2, col="red") +
  # improve viewing
  scale_y_continuous(breaks = seq(0,30,1)) + ylab("RND")
p2 = ggplot(plotInpFilter) +
  # TMRCA global
  geom_violin(aes(x=SpeciesPair, y=RND), alpha=0.1) +
  # TMRCA OF INVERSIONS
  geom_point(aes(x="bur-era", y=median(chr6$RND_bur_era)), alpha=0.6, pch=4, col="blue") +
  geom_point(aes(x="sar-num", y=median(chr6$RND_sar_num)), alpha=0.6, pch=4, col="blue") +
  geom_point(aes(x="bur-mel", y=median(chr6$RND_bur_mel)), alpha=0.6, pch=4, col="blue") +
  geom_point(aes(x="era-sar", y=median(chr6$RND_era_sar)), alpha=0.6, pch=4, col="blue") +
  geom_point(aes(x="mel-num", y=median(chr6$RND_mel_num)), alpha=0.6, pch=4, col="blue") +
  # TMRCA OF INVERSIONS - LINES
  geom_hline(aes(yintercept=median(chr6$RND_bur_mel)), alpha=0.6, lty=2, col="blue") +
  geom_hline(aes(yintercept=median(chr6$RND_era_sar)), alpha=0.6, lty=2, col="blue") +
  # improve viewing
  scale_y_continuous(breaks = seq(0,30,1)) + ylab("RND")
p3 = ggplot(plotInpFilter) +
  # TMRCA global
  geom_violin(aes(x=SpeciesPair, y=RND), alpha=0.1) +
  # TMRCA OF INVERSIONS
  geom_point(aes(x="bur-era", y=median(chr13$RND_bur_era)), alpha=0.6, pch=4, col="purple") +
  geom_point(aes(x="sar-num", y=median(chr13$RND_sar_num)), alpha=0.6, pch=4, col="purple") +
  geom_point(aes(x="bur-mel", y=median(chr13$RND_bur_mel)), alpha=0.6, pch=4, col="purple") +
  geom_point(aes(x="era-sar", y=median(chr13$RND_era_sar)), alpha=0.6, pch=4, col="purple") +
  geom_point(aes(x="mel-num", y=median(chr13$RND_mel_num)), alpha=0.6, pch=4, col="purple") +
  # TMRCA OF INVERSIONS - LINES
  geom_hline(aes(yintercept=median(chr13$RND_bur_era)), alpha=0.6, lty=2, col="purple") +
  geom_hline(aes(yintercept=median(chr13$RND_bur_mel)), alpha=0.6, lty=2, col="purple") +
  # improve viewing
  scale_y_continuous(breaks = seq(0,30,1)) + ylab("RND")
p4 = ggplot(plotInpFilter) +
  # TMRCA global
  geom_violin(aes(x=SpeciesPair, y=RND), alpha=0.1) +
  # TMRCA OF INVERSIONS
  geom_point(aes(x="bur-era", y=median(chr21$RND_bur_era)), alpha=0.6, pch=4, col="orange") +
  geom_point(aes(x="sar-num", y=median(chr21$RND_sar_num)), alpha=0.6, pch=4, col="orange") +
  geom_point(aes(x="bur-mel", y=median(chr21$RND_bur_mel)), alpha=0.6, pch=4, col="orange") +
  geom_point(aes(x="era-sar", y=median(chr21$RND_era_sar)), alpha=0.6, pch=4, col="orange") +
  geom_point(aes(x="mel-num", y=median(chr21$RND_mel_num)), alpha=0.6, pch=4, col="orange") +
  # TMRCA OF INVERSIONS - LINES
  geom_hline(aes(yintercept=median(chr21$RND_mel_num)), alpha=0.6, lty=2, col="orange") +
  # improve viewing
  scale_y_continuous(breaks = seq(0,30,1)) + ylab("RND")
g = grid.arrange(p1,p2,p3,p4, nrow=2)


#####





#################### TMRCA 95% intervals #########################
bur_era_025 = as.numeric(quantile(inpFilter$RND_bur_era, 0.05))
sar_num_025 = as.numeric(quantile(inpFilter$RND_sar_num, 0.05))
bur_mel_025 = as.numeric(quantile(inpFilter$RND_bur_mel, 0.05))
era_sar_025 = as.numeric(quantile(inpFilter$RND_era_sar, 0.05))
mel_num_025 = as.numeric(quantile(inpFilter$RND_mel_num, 0.05))

bur_era_975 = as.numeric(quantile(inpFilter$RND_bur_era, 0.95))
sar_num_975 = as.numeric(quantile(inpFilter$RND_sar_num, 0.95))
bur_mel_975 = as.numeric(quantile(inpFilter$RND_bur_mel, 0.95))
era_sar_975 = as.numeric(quantile(inpFilter$RND_era_sar, 0.95))
mel_num_975 = as.numeric(quantile(inpFilter$RND_mel_num, 0.95))

ggplot() +
  # INVERSIONS TMRCAs
  geom_boxplot(aes(x=0, y=10*chr2$RND_bur_mel ), alpha=0.5, fill="black", col="red") +
  geom_boxplot(aes(x=0, y=10*chr2$RND_era_sar ), alpha=0.5, fill="black", col="black") +
  # geom_boxplot(aes(x=1, y=10*chr6$RND_bur_mel ), alpha=0.5, fill="black", col="black") +
  # geom_boxplot(aes(x=2, y=10*chr6$RND_era_sar ), alpha=0.5, fill="black", col="black") +
  # geom_boxplot(aes(x=3, y=10*chr13$RND_bur_mel), alpha=0.5, fill="black", col="black") +
  # geom_boxplot(aes(x=4, y=10*chr13$RND_bur_era), alpha=0.5, fill="black", col="black") +
  # geom_boxplot(aes(x=5, y=10*chr21$RND_mel_num), alpha=0.5, fill="black", col="black") +
  # TMRCA 95% intervals
  annotate(geom = "rect", xmin=-1, xmax=6, ymin=10*bur_era_025, ymax=10*bur_era_975, fill="black", alpha=0.1) +
  annotate(geom = "rect", xmin=-1, xmax=6, ymin=10*sar_num_025, ymax=10*sar_num_975, fill="yellow", alpha=0.1) +
  annotate(geom = "rect", xmin=-1, xmax=6, ymin=10*bur_mel_025, ymax=10*bur_mel_975, fill="orange", alpha=0.1) +
  annotate(geom = "rect", xmin=-1, xmax=6, ymin=10*era_sar_025, ymax=10*era_sar_975, fill="red", alpha=0.1) +
  annotate(geom = "rect", xmin=-1, xmax=6, ymin=10*mel_num_025, ymax=10*mel_num_975, fill="blue", alpha=0.1) +
  # TMRCA 95% limits
  geom_hline(yintercept = c(bur_era_025,bur_era_975)*10, col="black") +
  geom_hline(yintercept = c(sar_num_025,sar_num_975)*10, col="yellow") +
  geom_hline(yintercept = c(bur_mel_025,bur_mel_975)*10, col="orange") +
  geom_hline(yintercept = c(era_sar_025,era_sar_975)*10, col="red") +
  geom_hline(yintercept = c(mel_num_025,mel_num_975)*10, col="blue")
  




ggplot(inpFilter) +
  # TMRCA global
  # between clades
  geom_boxplot(aes(x=1, y=RND_bur_era*10), fill="black", alpha=0.2) +
  geom_boxplot(aes(x=2, y=RND_sar_num*10), fill="purple", alpha=0.2) +
  geom_boxplot(aes(x=3, y=RND_bur_mel*10), fill="orange", alpha=0.2) +
  # within clades
  geom_boxplot(aes(x=4, y=RND_era_sar*10), fill="red", alpha=0.2) +
  geom_boxplot(aes(x=5, y=RND_mel_num*10), fill="blue", alpha=0.2) +
  # TMRCA of inversions
  # geom_hline(yintercept = chr21$RND_mel_num*10, col="green") +
  # geom_hline(yintercept = chr13$RND_bur_era*10, col="purple") +
  # geom_hline(yintercept = chr13$RND_bur_mel*10, col="purple") +
  geom_hline(yintercept = chr6$RND_bur_mel*10, col="red") +
  geom_hline(yintercept = chr6$RND_era_sar*10, col="pink") +
  # improve viewing
  scale_y_continuous(breaks = seq(0,30,1))



######################### Plot Data #########################
# plot
head(inpFilter)
ggplot(inpFilter) + 
  # geom_hline(yintercept = c(zscore.low.bur,zscore.hig.bur), col="red", alpha=0.5) +
  geom_line(aes(x=dummypos, y=RND_mel_num), col="blue", alpha=0.5) +
  geom_line(aes(x=dummypos, y=RND_era_sar), col="red", alpha=0.5) +
  geom_line(aes(x=dummypos, y=RND_bur_dor), col="orange", alpha=0.5) +
  scale_x_continuous(breaks = chromLabels$midpos, labels = seq(1,21,1)) +
  geom_vline(xintercept = c(min(chr2$dummypos), max(chr2$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr6$dummypos), max(chr6$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr13$dummypos), max(chr13$dummypos)), alpha=0.2, lty=2) +
  geom_vline(xintercept = c(min(chr21$dummypos), max(chr21$dummypos)), alpha=0.2, lty=2) +
  xlab("Chromosome") + ylab ("Relative Divergence") + ylim(0.0,2.2) + ggtitle("C") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))




# plot TRMCA of specific region vs global TMRCA distribution
# mel/num
p1 = ggplot(inpFilter) +
  geom_density(aes(RND_mel_num*10), fill="blue", alpha=0.2) +
  geom_vline(xintercept = median(chr2$RND_mel_num)*10, lty=2, col="black") +
  geom_vline(xintercept = median(chr6$RND_mel_num)*10, lty=2, col="red") +
  geom_vline(xintercept = median(chr13$RND_mel_num)*10, lty=2, col="purple") +
  geom_vline(xintercept = median(chr21$RND_mel_num)*10, lty=2, col="green") +
  ylim(0,1)
# era/sar
p2 = ggplot(inpFilter) +
  geom_density(aes(RND_era_sar*10), fill="red", alpha=0.2) +
  geom_vline(xintercept = median(chr2$RND_era_sar)*10, lty=2, col="black") +
  geom_vline(xintercept = median(chr6$RND_era_sar)*10, lty=2, col="red") +
  geom_vline(xintercept = median(chr13$RND_era_sar)*10, lty=2, col="purple") +
  geom_vline(xintercept = median(chr21$RND_era_sar)*10, lty=2, col="green") +
  ylim(0,1)
# bur/dor
p3 = ggplot(inpFilter) +
  geom_density(aes(RND_bur_dor*10), fill="orange", alpha=0.2) +
  geom_vline(xintercept = median(chr2$RND_bur_dor)*10, lty=2, col="black") +
  geom_vline(xintercept = median(chr6$RND_bur_dor)*10, lty=2, col="red") +
  geom_vline(xintercept = median(chr13$RND_bur_dor)*10, lty=2, col="purple") +
  geom_vline(xintercept = median(chr21$RND_bur_dor)*10, lty=2, col="green") +
  ylim(0,1)
# bur/mel
p4 = ggplot(inpFilter) +
  geom_density(aes(RND_bur_mel*10), fill="purple", alpha=0.2) +
  geom_vline(xintercept = median(chr2$RND_bur_mel)*10, lty=2, col="black") +
  geom_vline(xintercept = median(chr6$RND_bur_mel)*10, lty=2, col="red") +
  geom_vline(xintercept = median(chr13$RND_bur_mel)*10, lty=2, col="purple") +
  geom_vline(xintercept = median(chr21$RND_bur_mel)*10, lty=2, col="green") +
  ylim(0,1)
# bur/era
p5 = ggplot(inpFilter) +
  geom_density(aes(RND_bur_era*10), fill="pink", alpha=0.2) +
  geom_vline(xintercept = median(chr2$RND_bur_era)*10, lty=2, col="black") +
  geom_vline(xintercept = median(chr6$RND_bur_era)*10, lty=2, col="red") +
  geom_vline(xintercept = median(chr13$RND_bur_era)*10, lty=2, col="purple") +
  geom_vline(xintercept = median(chr21$RND_bur_era)*10, lty=2, col="green") +
  ylim(0,1)

grid.arrange(p1,p2,p3,p4,p5, nrow=1)





"era-mel"
"sar-num"
"bur-era"
"bur-mel"
"era-sar"
"mel-num"

# plot specific chromosome
inpFilterChrom = subset(inpFilter, scaffold == unique(chr21$scaffold))

ggplot(inpFilterChrom) + 
  # geom_hline(yintercept = c(zscore.low.bur,zscore.hig.bur), col="red", alpha=0.5) +
  geom_line(aes(x=dummypos, y=RND_mel_num), col="blue") +
  geom_line(aes(x=dummypos, y=RND_era_sar), col="red") +
  geom_line(aes(x=dummypos, y=RND_bur_dor), col="orange") +
  # geom_line(aes(x=dummypos, y=RND_bur_era), col="darkgreen") +
  # geom_line(aes(x=dummypos, y=RND_bur_mel), col="purple") +
  #
  geom_vline(xintercept = c(min(chr21$dummypos), max(chr21$dummypos)), alpha=0.2, lty=2) +
  xlab("Chromosome") + ylab ("Relative Divergence") + ylim(0.0,2.2) + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))

