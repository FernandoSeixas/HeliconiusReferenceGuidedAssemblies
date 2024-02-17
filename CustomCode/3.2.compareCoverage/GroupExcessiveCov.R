## load loibraries
require(ggplot2)
require(reshape2)
require(tidyr)
require(viridis)
require(gridExtra)

## variable
melClade = c("HEL_MEL","HEL_CYD","HEL_TIM","HEL_BES","HEL_NUM","HEL_HEC","HEL_ELE","HEL_PAR")
midClade = c("HEL_BUR","LAP_DOR")
eraClade = c("HEL_ERA","HEL_HIM_fat","HEL_HIM","HEL_SIA","HEL_TEL","HEL_DEM","HEL_SAR")

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

RepeatRegions = data.frame(
  chr = c("Hmel202001o","Hmel204001o","Hmel208001o","Hmel209001o","Hmel220003o","Hmel221001o"),
  sta = c(4075000,5650000,3300000,5125000,1050000,5825000),
  end = c(4125000,5875000,3475000,5450000,1075000,5850000),
  code = c("A","B","C","D","E","F"),
  plotNumber = c(1,2,3,4,5,6)
  )

# plot relative coverage for specific chromosome
sppTrio = c("HEL_HEC","HEL_ELE","HEL_PAR")
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


# heatmap
ch = bigchroms[18]
melCoverageSub = subset(melCoverage, scaffold == ch)
mmm = melCoverageSub[,c(1,2,3,6,7)]
ooo = dcast(mmm, scaffold + sta + end ~ reference)
ooo100k = data.frame()
for (rw in seq(1,nrow(ooo),4)) {
  s1 = ooo[rw:(rw+3),]
  ch = s1$scaffold[1]
  st = s1$sta[1]
  en = s1$end[4]
  s2 = colMeans(s1[,4:ncol(ooo)])
  df = data.frame(scaffold = ch, sta = st, end = en)
  ooo100k = rbind (ooo100k, cbind(df, t(s2)))
}
mat = t(as.matrix(ooo100k[,4:ncol(ooo100k)]))
heatmap(mat, Rowv = NA, Colv = NA, scale = "none", col = viridis(n = 256, direction = 1), main = ch)
# 
# ch = bigchroms[3]
# melCoverageSub = subset(melCoverage, scaffold == ch)
# subset(melCoverageSub, sta >= 375000 & end <= 400000 )
# # subset(melCoverageSub, relcov >= 5 & reference == "HEL_PAR")
# subset(melCoverageSub, sta >= 5000000 & end <= 6000000)

ch = unique(melCoverage$scaffold)[19]
melCoverageSub = subset(melCoverage, scaffold == ch )
xmaxim = max(melCoverageSub$end)/1000000
ggplot(melCoverageSub) + 
  geom_line(aes(x=(sta+end)/2000000, y=ypos+(relcov/max(melCoverageSub$relcov)), group=reference)) +
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


#################### regions of exceptional coverage in hec/ele/par ####################

# thresholds
minHig = 5
maxLow = 2

# simplify table
melCoverage$coordinate = paste0(melCoverage$scaffold,":",melCoverage$sta,":",melCoverage$end)
melCoverageTab = melCoverage[,c(13,6,7)]
melCoverageTab = dcast(melCoverageTab, coordinate ~ reference)
head(melCoverageTab)
melCoverageTab = melCoverageTab %>% separate(coordinate, c("scaffold","sta","end"), ":")

# transform to binary pattern
species = unique(melCoverage$reference)
melCoverageBin = melCoverageTab
for (cl in 4:ncol(melCoverageBin)) {
  melCoverageBin[,cl] = ifelse(melCoverageBin[,cl] > minHig, 2, melCoverageBin[,cl])
  melCoverageBin[,cl] = ifelse(melCoverageBin[,cl] < maxLow, 0, melCoverageBin[,cl])
  melCoverageBin[,cl] = ifelse(melCoverageBin[,cl] > maxLow & melCoverageBin[,cl] < minHig, 1, melCoverageBin[,cl])
}
melCoverageBin$sta = as.numeric(melCoverageBin$sta)
melCoverageBin$end = as.numeric(melCoverageBin$end)

## get patterns clade specific patterns
source("/n/home12/fseixas/code/heliconius_seixas/0.general/Rfunctions.R")

# sum binary values
melCoverageBin$out = (melCoverageBin$HEL_BES+melCoverageBin$HEL_BUR+melCoverageBin$HEL_CYD+melCoverageBin$HEL_DEM+melCoverageBin$HEL_ERA+melCoverageBin$HEL_HIM+melCoverageBin$HEL_HIM_fat+melCoverageBin$HEL_MEL+melCoverageBin$HEL_NUM+melCoverageBin$HEL_SAR+melCoverageBin$HEL_SIA+melCoverageBin$HEL_TEL+melCoverageBin$HEL_TIM)/2
melCoverageBin$ing = (melCoverageBin$HEL_ELE+melCoverageBin$HEL_PAR+melCoverageBin$HEL_HEC)/2
melCoverageBin$outmel = 
  (melCoverageBin$HEL_MEL+melCoverageBin$HEL_CYD+melCoverageBin$HEL_TIM+
     melCoverageBin$HEL_BES+melCoverageBin$HEL_NUM)/2
melCoverageBin$melClade = 
  (melCoverageBin$HEL_MEL+melCoverageBin$HEL_CYD+melCoverageBin$HEL_TIM+
  melCoverageBin$HEL_BES+melCoverageBin$HEL_NUM+
  melCoverageBin$HEL_HEC+melCoverageBin$HEL_ELE+melCoverageBin$HEL_PAR)/2
melCoverageBin$eraClade = 
  (melCoverageBin$HEL_ERA+melCoverageBin$HEL_HIM_fat+melCoverageBin$HEL_HIM+
   melCoverageBin$HEL_SIA+melCoverageBin$HEL_TEL+
   melCoverageBin$HEL_DEM+melCoverageBin$HEL_SAR)/2
  
subset(melCoverageBin, scaffold == bigchroms[19] & eraClade > 5)

mmm = subset(melCoverageBin, eraClade > 3 & melClade == 0)
mmmMerge = mergeRegions(mmm, 25000)
subset(mmmMerge, len > 25001)

ooo = subset(melCoverageBin, scaffold == bigchroms[19])


overlap = data.frame(matrix(nrow = 17, ncol = 17, data = 0))
for (x in 4:20) {
  for (y in 4:20) {
    n_A = nrow(melCoverageBin[melCoverageBin[,x] >= 1,])
    n_B = nrow(melCoverageBin[melCoverageBin[,y] > 1,])
    n_C = nrow(melCoverageBin)
    n_A_B = nrow(melCoverageBin[melCoverageBin[,x] > 1 & melCoverageBin[,y] > 1,])
    p = phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE) 
    # overlap[x-3,y-3] = n_A_B
    overlap[x-3,y-3] = p
  }
}
rownames(overlap) = names(melCoverageBin[4:20])
colnames(overlap) = names(melCoverageBin[4:20])
heatmap(as.matrix(overlap), Rowv = NA, Colv = NA, scale = "none", col = viridis(n = 256, direction = 1))

overlapLong = melt(overlap)
overlapLong$spp2 = rep(unique(overlapLong$variable), 17)
overlapLong = overlapLong[,c(1,3,2)]
names(overlapLong) = c("spp1","spp2","pvalue")



ggplot(overlapLong,aes(x=spp1,y=spp2,fill=pvalue))+
  geom_tile() +
  scale_fill_viridis()

subset(iii, pvalue > 0.01)
plot(density(iii$value))





# subset regions with outlier pattern
nrow(subset(melCoverageBin, ing >= 3 & out == 0)[,1:3])
HecEleParTrio = subset(melCoverageBin, ing >= 3 & out <= 0.0)[,1:3]
mergeRegions(HecEleParTrio, 25000)

nrow(subset(melCoverageBin, HEL_HEC >= 1 & HEL_ELE >= 1 & HEL_PAR >= 1 & ing >= 2.5 & outmel == 0)[,1:3])
HecEleParTrio2 = subset(melCoverageBin, HEL_HEC >= 1 & HEL_ELE >= 1 & HEL_PAR >= 1 & ing >= 2.5 & out == 0)[,1:3]
mergeRegions(HecEleParTrio2, 50000)


# return result
mrgAll$len = mrgAll$end - mrgAll$sta + 1

# melCoverageBin$sta = as.numeric(melCoverageBin$sta)
# melCoverageBin$end = as.numeric(melCoverageBin$end)
# ggplot(melCoverageBin) + 
#   # geom_line(aes(x=(sta+end)/2000000, y=melClade)) + 
#   geom_line(aes(x=(sta+end)/2000000, y=ing), col="blue") +
#   geom_line(aes(x=(sta+end)/2000000, y=melClade-ing), col="red") +
#   facet_wrap(~scaffold, scales="free_x", nrow=7)





subset(melCoverageBin, 
       HEL_MEL == 0 & HEL_CYD == 0 & HEL_TIM == 0 &
       HEL_PAR == 2 & HEL_ELE == 2 & HEL_HEC == 2 &
       HEL_BES == 0 & HEL_NUM == 0)[,1:3]


nrow(subset(melCoverageBin, HEL_PAR == 2 & HEL_ELE == 0 & HEL_HEC <= 1 & out == 0))
nrow(subset(melCoverageBin, HEL_PAR == 0 & HEL_ELE == 2 & HEL_HEC <= 1 & out == 0))
nrow(subset(melCoverageBin, ing >= 3 & HEL_HEC <= 1 & out == 0))
nrow(subset(melCoverageBin, HEL_PAR >= 2 & HEL_ELE >= 2 & HEL_HEC >= 2 & out == 0))

## plot along genomes
melCoverage$reference = factor(melCoverage$reference, levels = c(melClade,midClade,eraClade))
# all chromosomes
aaa = subset(melCoverage, reference %in% melClade)
ggplot(aaa) + 
  geom_line(aes(x=(sta+end)/2000000, y=relcov, col=reference)) + 
  facet_wrap(~scaffold, scales="free", nrow=7)
# specific chromosome
ch = unique(melCoverage$scaffold)[4]
aaa = subset(melCoverage, reference %in% melClade & scaffold == ch)
# ggplot(aaa) +
#   geom_area(aes(x=(sta+end)/2000000, y=relcov, fill=reference)) +
#   facet_wrap(~reference, scales="free", nrow=8) +
#   ggtitle(ch)
ggplot(aaa) +
  geom_area(aes(x=(sta+end)/2000000, y=relcov, fill=reference, colour=reference)) +
  facet_wrap(~reference, nrow=8) + xlim(5,6) + xlab("Chromosome Position (Mb)") +
  ggtitle(ch)


