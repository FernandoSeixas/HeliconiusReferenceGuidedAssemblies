## load loibraries
require(ggplot2)
require(reshape2)
require(tidyr)
require(viridis)
require(gridExtra)
require(stringr)


## variable
# melClade = c("HEL_MEL","HEL_CYD","HEL_TIM","HEL_BES","HEL_NUM","HEL_HEC","HEL_ELE","HEL_PAR")
# midClade = c("HEL_BUR","LAP_DOR")
# eraClade = c("HEL_ERA","HEL_HIM_fat","HEL_HIM","HEL_SIA","HEL_TEL","HEL_DEM","HEL_SAR")
melClade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
midClade = c("hbur","hdor")
eraClade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")


## get coverage statistics obtained by mapping to hmelv25 ==================================================
# dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudoreferences/4.2.compareMap/1.2.rmdup/hmelv25/"
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/3.1.compareMap/1.2.rmdup/hmelv25/"
files = list.files(path = dir, patter = ".cov.")
melCoverage = data.frame()
for (file in files) {
  inp = read.table(paste0(dir,file))
  names(inp) = c("scaffold","sta","end","sites","cov","reference")
  inp$relcov = inp$cov/median(inp$cov)
  inp$pos = seq(1,nrow(inp),1)
  inp$relpos = inp$pos/max(inp$pos)
  inp$set = "hmelv25"
  bigchroms = unique(subset(inp, end > 8000000)$scaffold)
  inpSub = subset(inp, scaffold %in% bigchroms)
  melCoverage = rbind(melCoverage, inpSub)
}
# rename species
melCoverage$reference = tolower(melCoverage$reference)
melCoverage$reference = str_replace(melCoverage$reference, "hel_","h")
melCoverage$reference = str_replace(melCoverage$reference, "lap_","h")
melCoverage$reference = str_replace(melCoverage$reference, "_", "")
# order species
melCoverage$reference = factor(melCoverage$reference, levels = c(melClade, midClade, eraClade))
# add clade info
melCoverage$clade = NA
melCoverage$clade = ifelse(melCoverage$reference %in% melClade, "melClade", melCoverage$clade)
melCoverage$clade = ifelse(melCoverage$reference %in% midClade, "midClade", melCoverage$clade)
melCoverage$clade = ifelse(melCoverage$reference %in% eraClade, "eraClade", melCoverage$clade)
# add ypos
melCoverage$ypos = -1
for (i in 1:length(levels(melCoverage$reference))) {
  spp = as.character(levels(melCoverage$reference)[i])
  melCoverage$ypos = ifelse(melCoverage$reference == spp, i, melCoverage$ypos)
}



######################### PLOT REGIONS OF INTEREST #########################
RepeatRegions = data.frame(
  chr = c("Hmel202001o","Hmel204001o","Hmel208001o","Hmel209001o","Hmel220003o","Hmel221001o"),
  sta = c(4075000,5650000,3300000,5125000,1050000,5825000),
  end = c(4125000,5875000,3475000,5450000,1075000,5850000),
  code = c("A","B","C","D","E","F"),
  plotNumber = c(1,2,3,4,5,6)
)
# plot relative coverage for specific chromosome
sppTrio = c("hhec","hele","hpar")
for (i in 1:6) {
  s1 = RepeatRegions[i,]
  # subset data for specific chromosome
  melCoverageSub = subset(melCoverage, scaffold == as.character(s1$chr) )
  melCoverageSub$trio = ifelse(melCoverageSub$reference %in% sppTrio, "Y", "N")
  # get limits of the repeat region and other specific info
  rlims = as.numeric(s1[1,2:3])
  ttlab = as.character(s1[1,4])
  pnumb = as.numeric(s1[1,5])
  xmaxim = floor(max(melCoverageSub$end)/1000000)
  # plot
  p = ggplot(melCoverageSub) + 
    geom_line(aes(x=(sta+end)/2000000, y=ypos+(relcov/max(melCoverageSub$relcov)), group=reference, col=trio)) +
    scale_x_continuous(breaks = seq(0,xmaxim,1)) +
    scale_colour_manual(values = c("darkgrey", "black")) +
    xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") +
    ggtitle(ttlab) +
    theme(
      plot.title = element_text(face="bold"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none" 
    )
  p = p + annotate("rect", xmin = (rlims[1]-25000)/1000000, xmax = (rlims[2]+25000)/1000000, ymin=0, ymax=Inf, fill="blue", alpha = 0.2)
  p
  # assign plot to variable
  assign(paste0("p", pnumb), p)
}
# combine plots
g = grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3)
ggsave(filename = paste0(dir,"../../RepeatRegions.png"), g, height = 10, width = 8, dpi = 150)




######################### REPEAT HOTSPOTS #########################
# specify unique coordinate codes
melCoverage$coordinate = paste0(melCoverage$scaffold,":",melCoverage$sta,":",melCoverage$end)

# transform relative coverage to binary
melCoverage$bin = ifelse(melCoverage$relcov > 2, 1, 0)

# count the number of species with excessice coverage (> 2x median coverage)
bins = aggregate(melCoverage$bin, list(melCoverage$coordinate), sum)
names(bins) = c("coordinate","count")
bins = bins %>% separate(coordinate, c("scaffold","sta","end"), ":")
bins$sta = as.numeric(bins$sta)
bins$end = as.numeric(bins$end)

# plot repeat hotspots
xmaxim = floor(max(melCoverage$end)/1000000)
a = ggplot(bins) +
  ## summary of outliers
  geom_tile(aes(x=(sta+end)/2000000, y=1, fill=count)) +
  scale_fill_viridis(direction = -1, name="Nb. of Individuals with\nRelative Coverage > 2") +
  facet_wrap(~scaffold, nrow=21, strip.position = "left") +
  scale_x_continuous(breaks = seq(0,xmaxim,1), expand = c(0, 0)) +
  scale_colour_manual(values = c("darkgrey", "black")) +
  xlab("Chromosome Position (Mb)") + ylab("Scaffold") +
  theme(
    plot.title = element_text(hjust = 0, vjust=0, face="bold"),
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    strip.text.y = element_text(angle=180),
    strip.background.y = element_rect(colour="white", fill="white"),
    strip.text.x = element_blank()
  ) +
  ggtitle("A")
ggsave(filename = paste0(dir, "RepeatHotspots.png"))

## plot repeat hotspots + relative coverage
ggplot(bins) +
  ## summary of outliers
  geom_tile(aes(x=(sta+end)/2000000, y=0, fill=count)) +
  scale_fill_viridis(direction = -1) +
  ## relative coverage per individual
  geom_line(data=melCoverage, aes(x=(sta+end)/2000000, y=ypos+(relcov/20), group=reference)) +
  scale_x_continuous(breaks = seq(0,xmaxim,1)) +
  scale_colour_manual(values = c("darkgrey", "black")) +
  xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") +
  ggtitle(ch) +
  theme(
    plot.title = element_text(hjust = 0, vjust=0, face="bold"),
    # legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  ) +
  facet_wrap(~scaffold, scales="free_x", ncol=7)
ggsave(filename = paste0(dir, "RelCoverage.png"), width = 210, height = 297, units = "mm")



######################### REPEAT HOTSPOTS RELATIVE POSITION #########################
# determine relative position in chromosomes
bins$RelPos = -1
for (ch in unique(bins$scaffold)) {
  s1 = subset(bins, scaffold == ch)
  lim = max(s1$end)
  bins$RelPos = ifelse(bins$scaffold == ch, bins$sta/lim, bins$RelPos)
}
bins$RelPos2 = abs(bins$RelPos - 0.5)
# plot  
b = ggplot(bins, aes(x=RelPos2, y=count)) + geom_smooth() + 
  xlab("Relative Distance to Chromosome Center") +
  ylab("Number of Species With\nExcessive Coverage at Window") +
  scale_y_continuous(breaks = seq(0,17,1)) +
  theme(
    plot.title = element_text(hjust = 0, vjust=0, face="bold"), panel.grid = element_blank()
  ) +
  ggtitle("C")



nrow(subset(bins, count > 3))/nrow(bins)*100

#################### regions of exceptional coverage in hec/ele/par ####################
# simplify table
melCoverageTab = melCoverage[,c(13,6,7)]
melCoverageTab = dcast(melCoverageTab, coordinate ~ reference)
melCoverageTab = melCoverageTab %>% separate(coordinate, c("scaffold","sta","end"), ":")

# transform to binary pattern
species = unique(melCoverage$reference)
melCoverageBin = melCoverageTab
for (cl in 4:ncol(melCoverageBin)) { melCoverageBin[,cl] = ifelse(melCoverageBin[,cl] > 2, 1, 0) }
melCoverageBin$sta = as.numeric(melCoverageBin$sta)
melCoverageBin$end = as.numeric(melCoverageBin$end)

# overlap of windows with excessive coverage (2x > median) 
overlap = data.frame(matrix(nrow = 17, ncol = 17, data = 0))
for (x in 4:20) {
  for (y in 4:20) {
    n_A = nrow(melCoverageBin[melCoverageBin[,x] == 1,])
    n_B = nrow(melCoverageBin[melCoverageBin[,y] == 1,])
    n_C = nrow(melCoverageBin)
    n_A_B = nrow(melCoverageBin[melCoverageBin[,x] == 1 & melCoverageBin[,y] == 1,])
    p = phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE) 
    overlap[x-3,y-3] = n_A_B
    # overlap[x-3,y-3] = p
  }
}
# 
names(overlap) = names(melCoverageBin)[4:20]
#
overlapLong = melt(overlap)
overlapLong$spp2 = rep(unique(overlapLong$variable), 17)
overlapLong = overlapLong[,c(1,3,2)]
#
names(overlapLong) = c("spp1","spp2","pvalue")
overlapLong$pvalue = round(overlapLong$pvalue, digits = 4)
# overlapLong$pvalue = ifelse(overlapLong$spp1 == overlapLong$spp2, NA, overlapLong$pvalue)
# plot
c = ggplot(overlapLong,aes(x=spp1,y=spp2,fill=pvalue))+
  geom_tile(col="white") +
  scale_fill_viridis(direction = -1, name="Nb. Overlaps") +
  xlab("") + ylab("") +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0, vjust=0, face="bold")
  ) +
  ggtitle("B")



grid.arrange(
  a,c,b,
  widths = c(4, 5, 1),
  layout_matrix = rbind(c(1, 1),
                        c(2, 3))
)
ggsave(filename = paste0(dir, "CNVsummary.png"), width = 210, height = 297, units = "mm")










# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# xdf = data.frame(Species = names(melCoverageBin)[4:20], Excessive = 0, PropExcessive = 0)
# for (x in 4:20) {
#   # xdf$Species = names(melCoverageBin)[x]
#   xdf$Excessive[(x-3)] = nrow(melCoverageBin[melCoverageBin[,x] >= 1,])
#   xdf$PropExcessive[(x-3)] = xdf$Excessive[(x-3)]/nrow(melCoverageBin)*100
# }
# 
# 
# 
# # subset regions with outlier pattern
# nrow(subset(melCoverageBin, ing >= 3 & out == 0)[,1:3])
# HecEleParTrio = subset(melCoverageBin, ing >= 3 & out <= 0.0)[,1:3]
# mergeRegions(HecEleParTrio, 25000)
# 
# nrow(subset(melCoverageBin, HEL_HEC >= 1 & HEL_ELE >= 1 & HEL_PAR >= 1 & ing >= 2.5 & outmel == 0)[,1:3])
# HecEleParTrio2 = subset(melCoverageBin, HEL_HEC >= 1 & HEL_ELE >= 1 & HEL_PAR >= 1 & ing >= 2.5 & out == 0)[,1:3]
# mergeRegions(HecEleParTrio2, 50000)
# 
# 
# # return result
# mrgAll$len = mrgAll$end - mrgAll$sta + 1
# 
# # melCoverageBin$sta = as.numeric(melCoverageBin$sta)
# # melCoverageBin$end = as.numeric(melCoverageBin$end)
# # ggplot(melCoverageBin) + 
# #   # geom_line(aes(x=(sta+end)/2000000, y=melClade)) + 
# #   geom_line(aes(x=(sta+end)/2000000, y=ing), col="blue") +
# #   geom_line(aes(x=(sta+end)/2000000, y=melClade-ing), col="red") +
# #   facet_wrap(~scaffold, scales="free_x", nrow=7)
# 
# 
# 
# 
# 
# subset(melCoverageBin, 
#        HEL_MEL == 0 & HEL_CYD == 0 & HEL_TIM == 0 &
#          HEL_PAR == 2 & HEL_ELE == 2 & HEL_HEC == 2 &
#          HEL_BES == 0 & HEL_NUM == 0)[,1:3]
# 
# 
# nrow(subset(melCoverageBin, HEL_PAR == 2 & HEL_ELE == 0 & HEL_HEC <= 1 & out == 0))
# nrow(subset(melCoverageBin, HEL_PAR == 0 & HEL_ELE == 2 & HEL_HEC <= 1 & out == 0))
# nrow(subset(melCoverageBin, ing >= 3 & HEL_HEC <= 1 & out == 0))
# nrow(subset(melCoverageBin, HEL_PAR >= 2 & HEL_ELE >= 2 & HEL_HEC >= 2 & out == 0))
# 
# ## plot along genomes
# melCoverage$reference = factor(melCoverage$reference, levels = c(melClade,midClade,eraClade))
# # all chromosomes
# aaa = subset(melCoverage, reference %in% melClade)
# ggplot(aaa) + 
#   geom_line(aes(x=(sta+end)/2000000, y=relcov, col=reference)) + 
#   facet_wrap(~scaffold, scales="free", nrow=7)
# # specific chromosome
# ch = unique(melCoverage$scaffold)[4]
# aaa = subset(melCoverage, reference %in% melClade & scaffold == ch)
# # ggplot(aaa) +
# #   geom_area(aes(x=(sta+end)/2000000, y=relcov, fill=reference)) +
# #   facet_wrap(~reference, scales="free", nrow=8) +
# #   ggtitle(ch)
# ggplot(aaa) +
#   geom_area(aes(x=(sta+end)/2000000, y=relcov, fill=reference, colour=reference)) +
#   facet_wrap(~reference, nrow=8) + xlim(5,6) + xlab("Chromosome Position (Mb)") +
#   ggtitle(ch)
# 
# 
