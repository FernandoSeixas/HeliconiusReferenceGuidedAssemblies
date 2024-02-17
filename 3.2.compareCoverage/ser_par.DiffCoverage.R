## load modules
require(reshape2)
require(ggplot2)
require(ggrepel)

dir = "/n/holyscratch01/mallet_lab/fseixas/2.elepar/1.2.rmdup/"
lfiles = list.files(dir, "cov.w25s25.txt")

## QTL regions
region1=2553006
region2=9375339



######################### READ DATA #########################
## read data
CNV = data.frame()
for (file in lfiles) {
  # read data
  inp = read.table(paste0(dir,file))
  names(inp) = c("scaffold","sta","end","sites","cov","ind")
  # coverage
  median(inp$cov)
  inp$relCoverage = inp$cov / median(inp$cov)
  # add to CNV
  CNV = rbind(CNV, inp)
}
## add pop specifications
CNV$pop = substr(CNV$ind, 1,12)

## Sex Dosage ////////////////////
# determine sex
s1 = subset(CNV, scaffold == "Hmel221001o")
zcov = aggregate(s1$relCoverage, list(s1$ind), median)
names(zcov) = c("ind", "zcov")
# apply sex dosage to sex chromosome
females = subset(zcov, zcov < 0.75)
for (i in unique(females$ind)) {
  CNV$relCoverage = ifelse(CNV$ind == i & CNV$scaffold == "Hmel221001o", CNV$relCoverage*2, CNV$relCoverage)
}



######################### PLOT "INTERACTION" #########################
# subset regions of interest
MaxSet1 = subset(CNV, scaffold == "Hmel221001o" & sta >= 2300000 & end <= 3125000)
MaxSet2 = subset(CNV, scaffold == "Hmel221001o" & sta >= 11025000 & end <= 11500000)
ggplot(MaxSet1) + 
  geom_point(aes(x=sta/1000000, y=relCoverage, col=pop, group=ind)) + 
  geom_vline(xintercept = 2300000/1000000) +
  geom_vline(xintercept = 2625000/1000000) +
  geom_vline(xintercept = c(2300000,2350000,3075000,3125000)/1000000, lty=2, col="grey")
  
ggplot(MaxSet2) + 
  geom_point(aes(x=sta/1000000, y=relCoverage, col=pop, group=ind)) +
  geom_vline(xintercept = 11025000/1000000) +
  geom_vline(xintercept = 11275000/1000000) +
  geom_vline(xintercept = 11450000/1000000) +
  geom_vline(xintercept = c(11025000,11075000,11350000,11500000)/1000000, lty=2, col="grey") 

MaxSet1 = subset(CNV, scaffold == "Hmel221001o" & sta == 2325000)
MaxSet2 = subset(CNV, scaffold == "Hmel221001o" & sta == 11275000)
# combine datasets
sub1 = MaxSet1[,6:8]
sub2 = MaxSet2[,6:8]
dfComb = merge(sub1,sub2, by = c("ind","pop"))
# plot
ggplot(dfComb, aes(x=relCoverage.x, y=relCoverage.y, col=as.factor(pop), group=pop, label=ind)) +
  geom_point(size=3, alpha=0.6) +
  scale_color_manual(values = c("#7570b3","#1b9e77"), name="Class") +
  geom_text_repel() +
  # xlim(2,10) + ylim(0,20) +
  xlab(paste0("Relative Coverage\n(Hmel221001o:2325000-2350000)")) +
  ylab("Relative Coverage\n(Hmel221001o:11450000-11475000)")



#################### COMPARE POPULATIONS #########################
popComparison = dcast(CNV[,c(1:3,6,7)], scaffold + sta + end ~ ind)
popComparison$ser = (popComparison$Hpar.ser.per.001+popComparison$Hpar.ser.per.002+popComparison$Hpar.ser.per.003+popComparison$Hpar.ser.per.004)/4
popComparison$but = (popComparison$Hpar.yur.per.001+popComparison$Hpar.yur.per.002+popComparison$Hpar.yur.per.003+popComparison$Hpar.yur.per.004)/4
popComparison$diff = popComparison$ser - popComparison$but

popComparison$diff.Zscore = (popComparison$diff - mean(popComparison$diff)) / sd(popComparison$diff)


butExcess = subset(popComparison, diff.Zscore < -3)[,1:3]
serExcess = subset(popComparison, diff.Zscore >  3)[,1:3]

source("/n/home12/fseixas/code/heliconius_seixas/0.general/Rfunctions.R")
butExcessMerge = mergeRegions(butExcess, 25000)
serExcessMerge = mergeRegions(serExcess, 25000)
nrow(butExcessMerge)
nrow(serExcessMerge)

merge(butExcess, CNV)

subset(butExcessMerge, len > 25001)
subset(serExcessMerge, len > 25001)



s1 = subset(popComparison, scaffold == "Hmel221001o")
maxcn = max(s1$relCoverage)
ggplot(s1) +
  geom_point(aes(x=(sta+end)/2000000, y=ser-but)) +
  geom_vline(xintercept = c(2300000,2350000)/1000000) +
  geom_vline(xintercept = c(3075000,3125000)/1000000) +
  # geom_vline(xintercept = c(11025000,11075000)/1000000, alpha=0.5) +
  # geom_vline(xintercept = c(11350000,11500000)/1000000, alpha=0.5) +
  scale_x_continuous(breaks = seq(0,30,1), limits=c(1,4)) + xlab("Chromosome Position (Mb)") +
  # scale_x_continuous(breaks = seq(0,30,1), limits=c(10,12)) + xlab("Chromosome Position (Mb)") +
  ggtitle(paste0("Hmel221001o",2300000,"-",3125000))
  # ggtitle(paste0("Hmel221001o",11350000,"-",11500000))
  


##### plot specific chromosome
s1 = subset(CNV, scaffold == "Hmel221001o")
maxcn = max(s1$relCoverage)
ggplot(s1) +
  geom_step(aes(x=(sta+end)/2000000, y=relCoverage, group=ind, col=pop)) +
  # geom_vline(xintercept = c(2300000,2350000)/1000000, alpha=0.5) +
  # geom_vline(xintercept = c(3075000,3125000)/1000000, alpha=0.5) +
  geom_vline(xintercept = c(11025000,11075000)/1000000, alpha=0.5) +
  geom_vline(xintercept = c(11350000,11500000)/1000000, alpha=0.5) +
  facet_wrap(~ind, nrow = 4) +
  # scale_x_continuous(breaks = seq(0,30,1), limits=c(1,4)) + xlab("Chromosome Position (Mb)") +
  scale_x_continuous(breaks = seq(0,30,1), limits=c(10,12)) + xlab("Chromosome Position (Mb)") +
  ylim(0,15) +
  ggtitle(paste0("Hmel221001o",2300000,"-",3125000)) +
  ggtitle(paste0("Hmel221001o",11350000,"-",11500000))








######################### ////////////////// #########################
head(CNV)
mat = dcast(CNV[,c(1:3,6,7)], scaffold + sta + end ~ ind)
mat$ser = (mat$Hpar.ser.per.001+mat$Hpar.ser.per.002+mat$Hpar.ser.per.003+mat$Hpar.ser.per.004)/4
mat$but = (mat$Hpar.yur.per.001+mat$Hpar.yur.per.002+mat$Hpar.yur.per.003+mat$Hpar.yur.per.004)/4


matsub$ser = (matsub$Hpar.ser.per.001+matsub$Hpar.ser.per.002+matsub$Hpar.ser.per.003+matsub$Hpar.ser.per.004)/4
matsub$but = (matsub$Hpar.yur.per.001+matsub$Hpar.yur.per.002+matsub$Hpar.yur.per.003+matsub$Hpar.yur.per.004)/4
ggplot(matsub) +
  geom_line(aes(x=(sta+end)/2000000, y=ser-but), col="black", alpha=0.5) +
  # geom_line(aes(x=(sta+end)/2000000, y=ser), alpha=0.8, col="red") +
  # geom_line(aes(x=(sta+end)/2000000, y=but), alpha=0.8, col="blue") +
  scale_x_continuous(breaks = seq(0,30,1)) + xlab("Chromosome Position (Mb)")


ggplot(matsub) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.ser.per.001/maxcn)+00), col="red", alpha=0.5) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.ser.per.002/maxcn)+01), col="red", alpha=0.5) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.ser.per.003/maxcn)+02), col="red", alpha=0.5) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.ser.per.004/maxcn)+03), col="red", alpha=0.7) +
  #
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.yur.per.001/maxcn)+04), col="black", alpha=0.7) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.yur.per.002/maxcn)+05), col="black", alpha=0.7) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.yur.per.003/maxcn)+06), col="black", alpha=0.7) +
  geom_line(aes(x=(sta+end)/2000000, y=(Hpar.yur.per.004/maxcn)+07), col="black", alpha=0.7) +
  #
  geom_vline(xintercept = c(region1,region2)/1000000, lty=2) +
  scale_x_continuous(breaks = seq(0,30,1)) + xlab("Chromosome Position (Mb)")



