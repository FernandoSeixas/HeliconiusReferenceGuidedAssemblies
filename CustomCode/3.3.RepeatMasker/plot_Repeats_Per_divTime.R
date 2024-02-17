## load libraries
require(ggplot2)
require(reshape2)

##
spp = "hhec"
dir = paste0("/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/",spp,".bigScaffs.Heliconius/")
max = 50

if (spp == "hmel") { WGsize = 268409364; RRsize = 325000 }
if (spp == "hhec") { WGsize = 299928006; RRsize = 3500000 }
if (spp == "hele") { WGsize = 323857173; RRsize = 4700000 }
if (spp == "hpar") { WGsize = 318856149; RRsize = 3300000 }


# ## REPEAT CLASSES ==================================================
# ## Repeat classes - WG
# fil = paste0(spp, "-2-hmelv25.HM.100gap.bigScaffs.fasta.out.landscape.Div.Rclass.tab")
# inp = read.table(paste0(dir,fil), header=T, skip=1)
# names(inp) = c("Rclass", paste0("[",seq(0,max-1,1),"-",seq(1,max,1),"]"))
# inpLong = melt(inp)
# sum(inpLong$value)
# xa = aggregate(inpLong$value, list(inpLong$Rclass), sum); names(xa) = c("Rclass","Count.WG")
# 
# a = ggplot(inpLong, aes(x=variable, y=value, fill=Rclass, col=Rclass)) +
#   geom_bar(position="stack", stat="identity", alpha=0.8) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0)) +
#   xlab("Percent Divergence from Concensus")
#   
# 
# ## Repeat classes - RR
# fil = paste0(spp, ".chr09_RepeatRegion.out.landscape.Div.Rclass.tab")
# inp = read.table(paste0(dir,fil), header=T, skip=1)
# names(inp) = c("Rclass", paste0("[",seq(0,max-1,1),"-",seq(1,max,1),"]"))
# inpLong = melt(inp)
# sum(inpLong$value)
# 
# xb = aggregate(inpLong$value, list(inpLong$Rclass), sum); names(xb) = c("Rclass","Count.RR")
# 
# b = ggplot(inpLong, aes(x=variable, y=value, fill=Rclass, col=Rclass)) +
#   geom_bar(position="stack", stat="identity", alpha=0.8) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0)) +
#   xlab("Percent Divergence from Concensus")
# 
# require(gridExtra)
# grid.arrange(a,b, nrow=2)
# xc = merge(xa,xb)
# xc$Prop.WG = round(xc$Count.WG/sum(xc$Count.WG)*100,2)
# xc$Prop.RR = round(xc$Count.RR/sum(xc$Count.RR)*100,2)
# 



## REPEAT FAMILIES ==================================================

## Repeat classes - WG
fil = paste0(spp, "-2-hmelv25.HM.100gap.bigScaffs.fasta.out.landscape.Div.Rfam.tab")
inpWG = read.table(paste0(dir,fil), header=T, skip=1)
# names(inp) = c("Rclass","Rfam",paste0("[",seq(0,max-1,1),"-",seq(1,max,1),"]"))
names(inpWG) = c("Rclass","Rfam",seq(1,max,1))

inpLongWG = melt(inpWG)
sum(inpLongWG$value)
xa = aggregate(inpLongWG$value, list(inpLongWG$Rfam), sum)
names(xa) = c("Rfam","Count.WG")

a = ggplot(inpLongWG, aes(x=variable, y=value/WGsize*100, fill=Rfam, col=Rfam)) +
  geom_bar(position="stack", stat="identity", alpha=0.8) +
  # scale_fill_viridis_d() +
  # scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0)) +
  xlab("Percent Divergence from Concensus") +
  ylab("Percent of Genome") 

## Repeat classes - RR
fil = paste0(spp, ".chr09_RepeatRegion.out.landscape.Div.Rfam.tab")
inpRR = read.table(paste0(dir,fil), header=T, skip=1)
names(inpRR) = c("Rclass","Rfam",seq(1,max,1))
inpLongRR = melt(inpRR)
sum(inpLongRR$value)
xb = aggregate(inpLongRR$value, list(inpLongRR$Rfam), sum)
names(xb) = c("Rfam","Count.RR")

b = ggplot(inpLongRR, aes(x=variable, y=value/RRsize*100, fill=Rfam, col=Rfam)) +
  geom_bar(position="stack", stat="identity", alpha=0.8) +
  # scale_fill_viridis_d() +
  # scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0)) +
  xlab("Percent Divergence from Concensus") +
  ylab("Percent of Repeat Region") +
  theme(legend.position = "none")

require(gridExtra)
grid.arrange(a,b, nrow=2)


##
agesWG = aggregate(inpLongWG$value, list(inpLongWG$variable), sum)
agesRR = aggregate(inpLongRR$value, list(inpLongRR$variable), sum)
names(agesWG) = c("Age","CountWG")
names(agesRR) = c("Age","CountRR")

agesWG$PropWG = agesWG$CountWG/sum(agesWG$CountWG)*100 
agesRR$PropRR = agesRR$CountRR/sum(agesRR$CountRR)*100 

agesComb = merge(agesWG[,c(1,3)],agesRR[,c(1,3)], sort = F)
agesCombLong = melt(agesComb)

ggplot(agesCombLong) +
  geom_col(aes(x=Age, y=value, fill=variable, col=variable), stat="identity", position="dodge", alpha=0.5) +
  # geom_line(aes(x=Age, y=value, group=variable, col=variable)) +
  scale_fill_manual(values = c("red","blue")) +
  scale_color_manual(values = c("red","blue")) +
  xlab("Percent Divergence from Concensus") +
  ylab("Proportion (%)") +
  ggtitle(spp)


##
xc = merge(xa,xb, by.x = "Rfam")
xc$Prop.WG = round(xc$Count.WG/sum(xc$Count.WG)*100,2)
xc$Prop.RR = round(xc$Count.RR/sum(xc$Count.RR)*100,2)
xc$Diff = xc$Prop.RR - xc$Prop.WG

ggplot(xc) + 
  geom_col(aes(x=Rfam, y=Diff), alpha=0.6, fill=ifelse(xc$Diff < 0, "red", "blue")) +
  ggtitle(spp)


inpLongYoungWG = subset(inpLongWG, as.numeric(variable) <= 5)
inpLongYoungRR = subset(inpLongRR, as.numeric(variable) <= 5)
xaa = aggregate(inpLongYoungWG$value, list(inpLongYoungWG$Rfam), sum)
xbb = aggregate(inpLongYoungRR$value, list(inpLongYoungRR$Rfam), sum)
names(xaa) = c("Rfam","Count.WG")
names(xbb) = c("Rfam","Count.RR")

xcc = merge(xaa,xbb, by.x = "Rfam")
xcc$Prop.WG = round(xcc$Count.WG/sum(xcc$Count.WG)*100,2)
xcc$Prop.RR = round(xcc$Count.RR/sum(xcc$Count.RR)*100,2)
xcc$Diff = xcc$Prop.RR - xcc$Prop.WG
ggplot(xcc) + 
  geom_col(aes(x=Rfam, y=Diff), alpha=0.6, fill=ifelse(xc$Diff < 0, "red", "blue")) +
  ggtitle(spp)
