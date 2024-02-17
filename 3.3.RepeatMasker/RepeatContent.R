## load libraries
require(ggplot2)

## fixed variables
chr = "chr08"
dir = "/n/scratchlfs/mallet_lab/fseixas/1.pseudoreferences/4.3.RepeatMasker/"

## define files
reffile = paste0(dir,"hmelv25.",chr,".target.Heliconius/",list.files(paste0(dir,"hmelv25.",chr,".target.Heliconius/"), pattern = ".transformed.out"))
hecfile = paste0(dir,"hhec.",chr,".target.Heliconius/",list.files(paste0(dir,"hhec.",chr,".target.Heliconius/"), pattern = ".transformed.out"))
elefile = paste0(dir,"hele.",chr,".target.Heliconius/",list.files(paste0(dir,"hele.",chr,".target.Heliconius/"), pattern = ".transformed.out"))
parfile = paste0(dir,"hpar.",chr,".target.Heliconius/",list.files(paste0(dir,"hpar.",chr,".target.Heliconius/"), pattern = ".transformed.out"))
# read files
ref = read.table(file=reffile, header=F, fill = T, row.names=NULL); ref = ref[!is.na(ref$V15),1:15]
hec = read.table(file=hecfile, header=F, fill = T, row.names=NULL); hec = hec[!is.na(hec$V15),1:15]
ele = read.table(file=elefile, header=F, fill = T, row.names=NULL); ele = ele[!is.na(ele$V15),1:15]
par = read.table(file=parfile, header=F, fill = T, row.names=NULL); par = par[!is.na(par$V15),1:15]
# header name
nms = c("bitScore","percDiv","percDel","PercIns","Query","Qsta","Qend","Qleft","symbol","matchRepeat","repeatFamily","Psta","Pend","repeatLeft","ID","star")
names(ref) = nms[1:(length(nms)-1)]
names(hec) = nms[1:(length(nms)-1)]
names(ele) = nms[1:(length(nms)-1)]
names(par) = nms[1:(length(nms)-1)]
#
ref$count = 1
hec$count = 1
ele$count = 1
par$count = 1
### count each repeat Family ==================================================
refCount = aggregate(ref$count, list(ref$repeatFamily), sum); names(refCount) = c("repeatFamily","Count"); refCount$relCount = refCount$Count/sum(refCount$Count); refCount$spp = "hmelv25"
hecCount = aggregate(hec$count, list(hec$repeatFamily), sum); names(hecCount) = c("repeatFamily","Count"); hecCount$relCount = hecCount$Count/sum(hecCount$Count); hecCount$spp = "hhec"
eleCount = aggregate(ele$count, list(ele$repeatFamily), sum); names(eleCount) = c("repeatFamily","Count"); eleCount$relCount = eleCount$Count/sum(eleCount$Count); eleCount$spp = "hele"
parCount = aggregate(par$count, list(par$repeatFamily), sum); names(parCount) = c("repeatFamily","Count"); parCount$relCount = parCount$Count/sum(parCount$Count); parCount$spp = "hpar"
# plot
cmbCount = rbind(refCount, hecCount, eleCount, parCount)
cmbCount$spp = factor(cmbCount$spp, levels = c("hmelv25","hhec","hele","hpar"))
lvls = refCount$repeatFamily[order(refCount$relCount, decreasing=T)]
cmbCount$repeatFamily = factor(cmbCount$repeatFamily, levels=lvls)

ggplot(cmbCount, aes(x=repeatFamily, y=Count, fill=spp, col=spp)) + 
  geom_col(position = "dodge", col="black", alpha=0.5) + 
  scale_fill_manual(values = c("black","orange","blue","red")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  facet_wrap(~spp, nrow=4) +
  xlab("Repeat Class/Family") + ylab("Count")

ggplot(cmbCount, aes(x=repeatFamily, y=relCount*100, fill=spp, col=spp)) + 
  geom_col(position = "dodge", col="black", alpha=0.5) + 
  scale_fill_manual(values = c("black","orange","blue","red")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  facet_wrap(~spp, nrow=4) +
  xlab("Repeat Class/Family") + ylab("Proportion (%)") + ggtitle(chr)

### Length each repeat Family ==================================================
if (chr == "chr08") {refTarget=175000; hecTarget=2100000; eleTarget=1900000; parTarget=2000000}
if (chr == "chr09") {refTarget=325000; hecTarget=4500000; eleTarget=4700000; parTarget=3300000}
refLen = aggregate(ref$Qend-ref$Qsta, list(ref$repeatFamily), sum); names(refLen) = c("repeatFamily","Length"); refLen$relLength = refLen$Length/refTarget; refLen$spp = "hmelv25"
hecLen = aggregate(hec$Qend-hec$Qsta, list(hec$repeatFamily), sum); names(hecLen) = c("repeatFamily","Length"); hecLen$relLength = hecLen$Length/hecTarget; hecLen$spp = "hhec"
eleLen = aggregate(ele$Qend-ele$Qsta, list(ele$repeatFamily), sum); names(eleLen) = c("repeatFamily","Length"); eleLen$relLength = eleLen$Length/eleTarget; eleLen$spp = "hele"
parLen = aggregate(par$Qend-par$Qsta, list(par$repeatFamily), sum); names(parLen) = c("repeatFamily","Length"); parLen$relLength = parLen$Length/parTarget; parLen$spp = "hpar"
# plot
cmbLen = rbind(refLen, hecLen, eleLen, parLen)
cmbLen$spp = factor(cmbLen$spp, levels = c("hmelv25","hhec","hele","hpar"))
lvls = refLen$repeatFamily[order(refLen$relLength, decreasing=T)]
cmbLen$repeatFamily = factor(cmbLen$repeatFamily, levels=lvls)

ggplot(cmbLen, aes(x=repeatFamily, y=relLength*100, fill=spp, col=spp)) + 
  geom_col(position = "dodge", col="black", alpha=0.5) + 
  scale_fill_manual(values = c("black","orange","blue","red")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  facet_wrap(~spp, nrow=4) +
  xlab("Repeat Class/Family") + ylab("Proportion of Region (%)") + ggtitle(chr)


# ### count each specific repeat ==================================================
# refCount = aggregate(ref$count, list(ref$matchRepeat), sum); names(refCount) = c("matchRepeat","Count"); refCount$relCount = refCount$Count/sum(refCount$Count); refCount$spp = "hmelv25"
# hecCount = aggregate(hec$count, list(hec$matchRepeat), sum); names(hecCount) = c("matchRepeat","Count"); hecCount$relCount = hecCount$Count/sum(hecCount$Count); hecCount$spp = "hhec"
# eleCount = aggregate(ele$count, list(ele$matchRepeat), sum); names(eleCount) = c("matchRepeat","Count"); eleCount$relCount = eleCount$Count/sum(eleCount$Count); eleCount$spp = "hele"
# parCount = aggregate(par$count, list(par$matchRepeat), sum); names(parCount) = c("matchRepeat","Count"); parCount$relCount = parCount$Count/sum(parCount$Count); parCount$spp = "hpar"
# # plot
# cmbCount = rbind(refCount, hecCount, eleCount, parCount)
# cmbCount$spp = factor(cmbCount$spp, levels = c("hmelv25","hhec","hele","hpar"))
# cmbCount = subset(cmbCount, relCount >= 0.01)
# ggplot(cmbCount, aes(x=matchRepeat, y=Count, fill=spp, col=spp)) + 
#   geom_col(position = "dodge", col="black", alpha=0.5) + 
#   scale_fill_manual(values = c("black","orange","blue","red")) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
#   facet_wrap(~spp, nrow=4) +
#   xlab("Repeat Class/Family") + ylab("Count")
# ggplot(cmbCount, aes(x=matchRepeat, y=relCount*100, fill=spp, col=spp)) + 
#   geom_col(position = "dodge", col="black", alpha=0.5) + 
#   scale_fill_manual(values = c("black","orange","blue","red")) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
#   facet_wrap(~spp, nrow=4) +
#   xlab("Repeat Class/Family") + ylab("Proportion (%)")



### genome background ==================================================
## define files
autreffile = paste0(dir,"hmel.autosomes.Heliconius/",list.files(paste0(dir,"hmel.autosomes.Heliconius/"), pattern = ".transformed.out"))
authecfile = paste0(dir,"hhec.autosomes.Heliconius/",list.files(paste0(dir,"hhec.autosomes.Heliconius/"), pattern = ".transformed.out"))
autelefile = paste0(dir,"hele.autosomes.Heliconius/",list.files(paste0(dir,"hele.autosomes.Heliconius/"), pattern = ".transformed.out"))
autparfile = paste0(dir,"hpar.autosomes.Heliconius/",list.files(paste0(dir,"hpar.autosomes.Heliconius/"), pattern = ".transformed.out"))
# read files
autref = read.table(file=autreffile, header=F, fill = T, row.names=NULL); autref = autref[!is.na(autref$V15),1:15]
authec = read.table(file=authecfile, header=F, fill = T, row.names=NULL); authec = authec[!is.na(authec$V15),1:15]
autele = read.table(file=autelefile, header=F, fill = T, row.names=NULL); autele = autele[!is.na(autele$V15),1:15]
autpar = read.table(file=autparfile, header=F, fill = T, row.names=NULL); autpar = autpar[!is.na(autpar$V15),1:15]
# header name
nms = c("bitScore","percDiv","percDel","PercIns","Query","Qsta","Qend","Qleft","symbol","matchRepeat","repeatFamily","Psta","Pend","repeatLeft","ID","star")
names(autref) = nms[1:(length(nms)-1)]
names(authec) = nms[1:(length(nms)-1)]
names(autele) = nms[1:(length(nms)-1)]
names(autpar) = nms[1:(length(nms)-1)]
#
autref$count = 1
authec$count = 1
autele$count = 1
autpar$count = 1
### count each repeat Family ==================================================
autrefCount = aggregate(autref$count, list(autref$repeatFamily), sum); names(autrefCount) = c("repeatFamily","Count"); autrefCount$relCount = autrefCount$Count/sum(autrefCount$Count); autrefCount$spp = "hmelv25"
authecCount = aggregate(authec$count, list(authec$repeatFamily), sum); names(authecCount) = c("repeatFamily","Count"); authecCount$relCount = authecCount$Count/sum(authecCount$Count); authecCount$spp = "hhec"
auteleCount = aggregate(autele$count, list(autele$repeatFamily), sum); names(auteleCount) = c("repeatFamily","Count"); auteleCount$relCount = auteleCount$Count/sum(auteleCount$Count); auteleCount$spp = "hele"
autparCount = aggregate(autpar$count, list(autpar$repeatFamily), sum); names(autparCount) = c("repeatFamily","Count"); autparCount$relCount = autparCount$Count/sum(autparCount$Count); autparCount$spp = "hpar"
# plot
autcmbCount = rbind(autrefCount, authecCount, auteleCount, autparCount)
autcmbCount$spp = factor(autcmbCount$spp, levels = c("hmelv25","hhec","hele","hpar"))
lvls = autrefCount$repeatFamily[order(autrefCount$relCount, decreasing=T)]
autcmbCount$repeatFamily = factor(autcmbCount$repeatFamily, levels=lvls)

ggplot(autcmbCount, aes(x=repeatFamily, y=Count, fill=spp, col=spp)) + 
  geom_col(position = "dodge", col="black", alpha=0.5) + 
  scale_fill_manual(values = c("black","orange","blue","red")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  facet_wrap(~spp, nrow=4) +
  xlab("Repeat Class/Family") + ylab("Count")
ggplot(autcmbCount, aes(x=repeatFamily, y=relCount*100, fill=spp, col=spp)) + 
  geom_col(position = "dodge", col="black", alpha=0.5) + 
  scale_fill_manual(values = c("black","orange","blue","red")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  facet_wrap(~spp, nrow=4) +
  xlab("Repeat Class/Family") + ylab("Proportion (%)") + ggtitle(chr)


## Compare
MrgCount = merge(cmbCount, autcmbCount, by = c("repeatFamily","spp"))
MrgCount$diff = MrgCount$relCount.x-MrgCount$relCount.y
lvls = unique(MrgCount$repeatFamily[order(MrgCount$relCount.x-MrgCount$relCount.y, decreasing = T)])
MrgCount$repeatFamily = factor(MrgCount$repeatFamily, levels=lvls)

ggplot(MrgCount, aes(x=repeatFamily, y=diff*100, fill=spp, col=spp)) + 
  geom_col(position = "dodge", col="black", alpha=0.5) + 
  scale_fill_manual(values = c("black","orange","blue","red")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  facet_wrap(~spp, nrow=4) +
  xlab("Repeat Class/Family") + ylab("Proportion (%)") + ggtitle(chr)
