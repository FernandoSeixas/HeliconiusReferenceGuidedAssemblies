## load libraries
require(ggplot2)

## read data
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.StructuralVariation/covFluctuations/mel2trio/3.coverage/"
hecfil = paste0("HEL_MEL_b1.to.hhecMEDUSA.coverage.txt")
elefil = paste0("HEL_MEL_b1.to.heleMEDUSA.coverage.txt")
parfil = paste0("HEL_MEL_b1.to.hparMEDUSA.coverage.txt")
hec = read.table(paste0(dir,hecfil))
ele = read.table(paste0(dir,elefil))
par = read.table(paste0(dir,parfil))
names(hec) = c("scaffold","sta","end","reads","coverage","ind")
names(ele) = c("scaffold","sta","end","reads","coverage","ind")
names(par) = c("scaffold","sta","end","reads","coverage","ind")
hecmed = median(hec$coverage)
elemed = median(ele$coverage)
parmed = median(par$coverage)

## 
cnb = 9
subhec = subset(hec, scaffold == paste0("hhec000",cnb))
subele = subset(ele, scaffold == paste0("hele000",cnb))
subpar = subset(par, scaffold == paste0("hpar000",cnb))
ggplot(subhec, aes(x=(sta+end)/2000000, y=coverage/hecmed)) + geom_point(alpha=0.5, col="orange") + ylim(0,5) + scale_x_continuous(breaks = seq(0,20,1.0)) + xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle(paste0("hhec000",cnb))
ggplot(subele, aes(x=(sta+end)/2000000, y=coverage/elemed)) + geom_point(alpha=0.5, col="blue")   + ylim(0,5) + scale_x_continuous(breaks = seq(0,20,1.0)) + xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle(paste0("hele000",cnb))
ggplot(subpar, aes(x=(sta+end)/2000000, y=coverage/parmed)) + geom_point(alpha=0.5, col="red")    + ylim(0,5) + scale_x_continuous(breaks = seq(0,20,1.0)) + xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle(paste0("hpar000",cnb))
ggplot(subhec, aes(x=(sta+end)/2000000, y=coverage/hecmed)) + geom_line(alpha=0.5, col="orange") + ylim(0,4) + scale_x_continuous(limits=c(3.0,6.0), breaks = seq(0,20,0.2)) + xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle(paste0("hhec000",cnb))
ggplot(subele, aes(x=(sta+end)/2000000, y=coverage/elemed)) + geom_line(alpha=0.5, col="blue")   + ylim(0,4) + scale_x_continuous(limits=c(3.0,6.0), breaks = seq(0,20,0.2)) + xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle(paste0("hele000",cnb))
ggplot(subpar, aes(x=(sta+end)/2000000, y=coverage/parmed)) + geom_line(alpha=0.5, col="red")    + ylim(0,4) + scale_x_continuous(limits=c(3.0,6.0), breaks = seq(0,20,0.2)) + xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle(paste0("hpar000",cnb))

# hhec0008 borders: 3.2-5.3 Mb -- 2.1 Mb
# hele0008 borders: 3.5-5.4 Mb -- 1.9 Mb
# hpar0008 borders: 3.6-5.6 Mb -- 2.0 Mb

# hhec0009 borders: 5.5-9.0  Mb -- 4.5 Mb
# hele0009 borders: 6.5-11.2 Mb -- 4.7.Mb
# hpar0009 borders: 6.5-9.8  Mb -- 3.3 Mb
