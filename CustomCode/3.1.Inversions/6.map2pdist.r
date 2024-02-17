## load modules
require(ggplot2)
require(gridExtra)

## variables
chr="Hmel202001o"; bk1=2147755; bk2=2230606
chr="Hmel206001o"; bk1=12129405; bk2=12591271; bk3=12279100; bk4=12280171
chr="Hmel213001o"; bk1=16917575; bk2=17154514
chr="Hmel221001o"; bk1=8021528; bk2=8342578

## read data  ==================================================
ref = "hmelv25"
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/fdM/5.pdistance/"
fil = paste0(chr,".whole.outgroup.dist.csv")
pdistance = read.csv(paste0(dir,fil))

pdistanceFilter = subset(pdistance, sites >= 5000)



## plot 
# ggplot(pdistanceFilter) +
#   geom_line(aes(x=(start+end)/2000000, y=dxy_era_bur/dxy_era_mel), col="red", lty=1) +
#   geom_line(aes(x=(start+end)/2000000, y=dxy_era_dor/dxy_era_mel), col="red", lty=2) +
#   geom_line(aes(x=(start+end)/2000000, y=dxy_mel_bur/dxy_era_mel), col="blue", lty=1) +
#   geom_line(aes(x=(start+end)/2000000, y=dxy_mel_dor/dxy_era_mel), col="blue", lty=2) +
#   geom_vline(xintercept = c(bk1,bk2)/1000000) +
#   xlab("Chromosome Position (Mb)") + ylab("dxy") + ylim(0,2)

## RND
p1 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_era_bur/dxy_era_mel), col="#d95f02") +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,2) + 
  ggtitle("hera-hbur") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p2 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_era_dor/dxy_era_mel), col="#d95f02", alpha=0.8) +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,2) + 
  ggtitle("hera-hdor") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
p3 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_mel_bur/dxy_era_mel), col="#7570b3") +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,2) + 
  ggtitle("hmel-hbur") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p4 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_mel_dor/dxy_era_mel), col="#7570b3", alpha=0.8) +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,2) + 
  ggtitle("hmel-hdor") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
grid.arrange(p1, p3, p2, p4, nrow=2)  


## simple dxy 
p1 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_era_bur), col="#d95f02") +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("DXY") + ylim(0,0.1) + 
  ggtitle("hera-hbur") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p2 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_era_dor), col="#d95f02", alpha=0.8) +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("DXY") + ylim(0,0.1) + 
  ggtitle("hera-hdor") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p3 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_mel_bur), col="#7570b3") +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("DXY") + ylim(0,0.1) + 
  ggtitle("hmel-hbur") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p4 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_mel_dor), col="#7570b3", alpha=0.8) +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("DXY") + ylim(0,0.1) + 
  ggtitle("hmel-hdor") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
grid.arrange(p1, p3, p2, p4, nrow=2)  

#
p1 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_era_bur/dxy_bur_dor), col="#d95f02") +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,2) + 
  ggtitle("hera-hbur") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p2 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_era_dor/dxy_bur_dor), col="#d95f02", alpha=0.8) +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,2) + 
  ggtitle("hera-hdor") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p3 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_mel_bur/dxy_bur_dor), col="#7570b3") +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,1.5) + 
  ggtitle("hmel-hbur") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
p4 <- ggplot(pdistanceFilter) + 
  geom_line(aes(x=(start+end)/2000000, y=dxy_mel_dor/dxy_bur_dor), col="#7570b3", alpha=0.8) +
  geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
  xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0.5,1.5) + 
  ggtitle("hmel-hdor") + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
grid.arrange(p1, p3, p2, p4, nrow=2)  



quantile(pdistanceFilter$dxy_era_bur/pdistanceFilter$dxy_era_mel, 0.01)
mmm = subset(pdistanceFilter, end >= bk1 & start <= bk2)

min(mmm$dxy_era_bur/mmm$dxy_era_mel)



# ## RND
# p1 <- ggplot(pdistanceFilter) + 
#   geom_line(aes(x=(start+end)/2000000, y=dxy_era_bur/((dxy_bur_out+dxy_dor_out)/2)), col="#d95f02") +
#   geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
#   xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0,1) + 
#   ggtitle("hera-hbur") + theme_bw() + 
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
# p2 <- ggplot(pdistanceFilter) + 
#   geom_line(aes(x=(start+end)/2000000, y=dxy_era_dor/((dxy_bur_out+dxy_dor_out)/2)), col="#d95f02", alpha=0.8) +
#   geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
#   xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0,1) + 
#   ggtitle("hera-hdor") + theme_bw() + 
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
# p3 <- ggplot(pdistanceFilter) + 
#   geom_line(aes(x=(start+end)/2000000, y=dxy_mel_bur/((dxy_bur_out+dxy_dor_out)/2)), col="#7570b3") +
#   geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
#   xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0,1) + 
#   ggtitle("hmel-hbur") + theme_bw() + 
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
# p4 <- ggplot(pdistanceFilter) + 
#   geom_line(aes(x=(start+end)/2000000, y=dxy_mel_dor/((dxy_bur_out+dxy_dor_out)/2)), col="#7570b3", alpha=0.8) +
#   geom_vline(xintercept = c(bk1,bk2)/1000000, alpha=0.2) + 
#   xlab("Chromosome Position (Mb)") + ylab("RND") + ylim(0,1) + 
#   ggtitle("hmel-hdor") + theme_bw() + 
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12))
# grid.arrange(p1, p3, p2, p4, nrow=2)  
