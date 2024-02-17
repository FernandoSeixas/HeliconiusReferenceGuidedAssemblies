require(reshape2)
require(ggplot2)

## variables
ref="hmelv25"
spp="hele"
max = 50

refWGsize = 268409364; refRRsize = 325000
if (spp == "hhec") { sppWGsize = 299928006; sppRRsize = 3500000 }
if (spp == "hele") { sppWGsize = 323857173; sppRRsize = 4700000 }
if (spp == "hpar") { sppWGsize = 318856149; sppRRsize = 3300000 }

dir = paste0("/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/")

refLfiles = list.files(pattern = "Rfam.tab", paste0(dir,ref,".bigScaffs.Heliconius"))
sppLfiles = list.files(pattern = "Rfam.tab", paste0(dir,spp,".bigScaffs.Heliconius"))


##
refWG = read.table(paste0(dir,ref,".bigScaffs.Heliconius/",refLfiles[!grepl(pattern = "RepeatRegion", refLfiles)]), header=T, skip=1)
sppWG = read.table(paste0(dir,spp,".bigScaffs.Heliconius/",sppLfiles[!grepl(pattern = "RepeatRegion", sppLfiles)]), header=T, skip=1)
refRR = read.table(paste0(dir,ref,".bigScaffs.Heliconius/",refLfiles[grepl(pattern = "RepeatRegion", refLfiles)]), header=T, skip=1)
sppRR = read.table(paste0(dir,spp,".bigScaffs.Heliconius/",sppLfiles[grepl(pattern = "RepeatRegion", sppLfiles)]), header=T, skip=1)

names(refWG) = c("Rclass","Rfam",seq(1,max,1))
names(sppWG) = c("Rclass","Rfam",seq(1,max,1))
names(refRR) = c("Rclass","Rfam",seq(1,max,1))
names(sppRR) = c("Rclass","Rfam",seq(1,max,1))

refWGLong = melt(refWG)
sppWGLong = melt(sppWG)
refRRLong = melt(refRR)
sppRRLong = melt(sppRR)

refWGagg = aggregate(refWGLong$value, list(refWGLong$Rfam), sum)
sppWGagg = aggregate(sppWGLong$value, list(sppWGLong$Rfam), sum)
refRRagg = aggregate(refRRLong$value, list(refRRLong$Rfam), sum)
sppRRagg = aggregate(sppRRLong$value, list(sppRRLong$Rfam), sum)

sum(refWGagg[,2])/refWGsize*100
sum(sppWGagg[,2])/sppWGsize*100
sum(refRRagg[,2])/refRRsize*100
sum(sppRRagg[,2])/sppRRsize*100


## combined data.frame
rfams = unique(
  c(
    as.character(refWGagg$Group.1),
    as.character(sppWGagg$Group.1),
    as.character(refRRagg$Group.1),
    as.character(sppRRagg$Group.1))
)

rfamDF = data.frame(Rfam=rfams, refWGcount = 0, sppWGcount = 0, refRRcount = 0, sppRRcount = 0)
for (rw in 1:nrow(rfamDF)) {
  rf = as.character(rfamDF[rw,1])
  if ( length(refWGagg[refWGagg$Group.1 == rf,2]) >0 )  { rfamDF$refWGcount[rw] = refWGagg[refWGagg$Group.1 == rf,2] } 
  if ( length(sppWGagg[sppWGagg$Group.1 == rf,2]) >0 )  { rfamDF$sppWGcount[rw] = sppWGagg[sppWGagg$Group.1 == rf,2] } 
  if ( length(refRRagg[refRRagg$Group.1 == rf,2]) >0 )  { rfamDF$refRRcount[rw] = refRRagg[refRRagg$Group.1 == rf,2] } 
  if ( length(sppRRagg[sppRRagg$Group.1 == rf,2]) >0 )  { rfamDF$sppRRcount[rw] = sppRRagg[sppRRagg$Group.1 == rf,2] } 
  # if ( length(refWGagg[refWGagg$Group.1 == rf,2]) >0 )  { rfamDF$refWGcount[rw] = refWGagg[refWGagg$Group.1 == rf,2]/refWGsize*100 } 
  # if ( length(sppWGagg[sppWGagg$Group.1 == rf,2]) >0 )  { rfamDF$sppWGcount[rw] = sppWGagg[sppWGagg$Group.1 == rf,2]/sppWGsize*100 } 
  # if ( length(refRRagg[refRRagg$Group.1 == rf,2]) >0 )  { rfamDF$refRRcount[rw] = refRRagg[refRRagg$Group.1 == rf,2]/refRRsize*100 } 
  # if ( length(sppRRagg[sppRRagg$Group.1 == rf,2]) >0 )  { rfamDF$sppRRcount[rw] = sppRRagg[sppRRagg$Group.1 == rf,2]/sppRRsize*100 } 
}

rfamDF$WGdiff = rfamDF$sppWGcount - rfamDF$refWGcount
rfamDF$RRdiff = rfamDF$sppRRcount - rfamDF$refRRcount
rfamDFlong = melt(rfamDF[,c(1,6:7)])

ggplot(rfamDFlong) +
  geom_bar(aes(x=Rfam, y=value, fill=variable, group=variable), stat = "identity", position="dodge") +
  ggtitle(paste0(ref, "-vs-", spp))


# 
# rfamDFlong = melt(rfamDF[,1:5])
# ggplot(rfamDFlong) +
#   geom_bar(aes(x=Rfam, y=value, fill=variable, group=variable), stat = "identity", position="dodge")
# 
# 
# for (rw in 1:nrow(rfamDF)) {
#   df <- matrix(as.numeric(rfamDF[rw,2:5]), ncol = 2, dimnames = list(c("ref", "spp"), c("wg", "rr")))
#   fet = fisher.test(df, alternative = "greater")
#   rfamDF$fet[rw] = fet$p.value
# }
# 
# rfamDF$sign = ifelse(rfamDF$fet < 0.01, 1, 0)
# options(scipen = 1)
# rfamDF
