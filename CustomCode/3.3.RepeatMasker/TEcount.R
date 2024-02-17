## load modules
require(stringr)
require(tidyr)
require(ggplot2)
require(reshape2)


##### read data ==================================================
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/TEcounts/"

lfiles = list.files(path = dir, pattern = ".Rfam.")
lfiles = lfiles[!grepl(pattern = "RR02", lfiles)]
lfiles = lfiles[!grepl(pattern = "RR08", lfiles)]
# lfiles = lfiles[!grepl(pattern = "RR09", lfiles)]

df = data.frame() 
for (fl in lfiles) {
  inp = read.table(paste0(dir,fl))
  # inp$prop = inp$V2/sum(inp$V2)
  inp$set = fl
  df = rbind(df, inp)
}
names(df) = c("Rfam","Count","set")
# names(df) = c("Rfam","Count","Prop","set")

# define species and set [wg, rr] in different columns 
df$set = str_remove(df$set, ".Rfam")
df$set = str_remove(df$set, ".txt")
df = separate(data=df, col=set, into = c("Species","Region"), "[.]")
df$set = paste0(df$Species,".",df$Region)

## remove RR counts from WG
dfCast = dcast(df, Rfam ~ set, value.var = "Count")
dfCast[is.na(dfCast)] = 0
head(dfCast)

# remove RR counts
dfCast$hmelv25.WG = dfCast$hmelv25.WG - dfCast$hmelv25.RR09
dfCast$hhec.WG = dfCast$hhec.WG - dfCast$hhec.RR09
dfCast$hele.WG = dfCast$hele.WG - dfCast$hele.RR09
dfCast$hpar.WG = dfCast$hpar.WG - dfCast$hpar.RR09

# melt data 
df = melt(dfCast, id.vars = "Rfam")
names(df) = c("Rfam","set","Count")

# add proportion
df$Prop = 0
for (i in unique(df$set)) {
  total = sum(subset(df, set == i)$Count)
  df$Prop = ifelse(df$set == i, df$Count/total, df$Prop)
}

# define species and set [wg, rr] in different columns 
df$set = str_remove(df$set, ".Rfam")
df$set = str_remove(df$set, ".txt")
df = separate(data=df, col=set, into = c("Species","Region"), "[.]")
df$set = paste0(df$Species,".",df$Region)

## plot all Rfams per species+region
ggplot(df) +
  geom_bar(aes(x=Rfam, y=Prop, group=set, fill=set, col=set), stat="identity", position="dodge", alpha=0.6) +
  coord_flip()


############### Enrichment vs. melpomene ? ##################
################## Compare Proportions ######################
#// this tests if the relative number of repeats (s or g) in the repeat region in relation to the whole genome (S or G) is
#// larger in the target species (hhec,hele,hpar) than in the reference species (hmel)
#// (s/S)/(g/G)

# transform data
dfCast = dcast(df, Rfam ~ set, value.var = "Prop")
dfCast[is.na(dfCast)] = 0

# relative proportion of each Rfam in target species (RR/WG) versus hmelv25 (RR/WG)
dfCast$hhecDprop = (dfCast$hhec.RR/dfCast$hhec.WG)/(dfCast$hmelv25.RR/dfCast$hmelv25.WG)
dfCast$heleDprop = (dfCast$hele.RR/dfCast$hele.WG)/(dfCast$hmelv25.RR/dfCast$hmelv25.WG)
dfCast$hparDprop = (dfCast$hpar.RR/dfCast$hpar.WG)/(dfCast$hmelv25.RR/dfCast$hmelv25.WG)

# only finite values
dfCast = subset(dfCast, is.finite(dfCast[,10]) & is.finite(dfCast[,11]) & is.finite(dfCast[,12]))

# enriched in the RR of the trio
subset(dfCast,hhecDprop > 1 & hparDprop > 1 & hparDprop > 1 )

# #
# dfCast$hmelDiff = (dfCast$hmelv25.RR-dfCast$hmelv25.WG)*100
# dfCast$hhecDiff = (dfCast$hhec.RR-dfCast$hhec.WG)*100
# dfCast$heleDiff = (dfCast$hele.RR-dfCast$hele.WG)*100
# dfCast$hparDiff = (dfCast$hpar.RR-dfCast$hpar.WG)*100
# 
# dfCastLong = melt(dfCast[,c(1,13:16)])
# 
# # ggplot(subset(dfCastLong, variable %in% c("hmelDiff","hparDiff"))) +
# ggplot(dfCastLong) +
#   geom_bar(aes(x=Rfam, y=value, fill=variable), stat="identity", position="dodge") +
#   scale_fill_brewer(type = "qual", palette = 6) +
#   coord_flip() 
# 
# subset(dfCast, hmelDiff <= 1 & hhecDiff > 1 & heleDiff > 1 & hparDiff > 1)
# 
# subset(dfCast, 
#        hhecDprop > 1 & hparDprop > 1 & hparDprop > 1 &
#        hmelDiff <= 1 & hhecDiff > 1 & heleDiff > 1 & hparDiff > 1
#        )
# 
# 

############### Enrichment vs. melpomene ? ################
#################### Fisher exact test ####################

dfCast = dcast(df, Rfam ~ set, value.var = "Count")
dfCast[is.na(dfCast)] = 0

dfCast$fetHec = 0
dfCast$fetEle = 0
dfCast$fetPar = 0

for (rw in 1:nrow(dfCast)) {
  # rfam
  print(dfCast$Rfam[rw])
  # counts
  refRR = dfCast$hmelv25.RR[rw]; refWG = dfCast$hmelv25.WG[rw]
  hecRR = dfCast$hhec.RR[rw]; hecWG = dfCast$hhec.WG[rw]
  eleRR = dfCast$hele.RR[rw]; eleWG = dfCast$hele.WG[rw]
  parRR = dfCast$hpar.RR[rw]; parWG = dfCast$hpar.WG[rw]
  # matrices for FET
  mat_hec = matrix(c(refWG,hecWG,refRR,hecRR), ncol=2, dimnames = list(c("ref", "spp"), c("wg", "rr")))
  mat_ele = matrix(c(refWG,eleWG,refRR,eleRR), ncol=2, dimnames = list(c("ref", "spp"), c("wg", "rr")))
  mat_par = matrix(c(refWG,parWG,refRR,parRR), ncol=2, dimnames = list(c("ref", "spp"), c("wg", "rr")))
  # FET
  fet_hec = fisher.test(mat_hec, alternative = "greater")
  fet_ele = fisher.test(mat_ele, alternative = "greater")
  fet_par = fisher.test(mat_par, alternative = "greater")
  # add info to table
  dfCast$fetHec[rw] = fet_hec$p.value
  dfCast$fetEle[rw] = fet_ele$p.value
  dfCast$fetPar[rw] = fet_par$p.value
}


nrow(subset(dfCast, fetHec < 0.01 & fetEle < 0.01 & fetPar < 0.01))
subset(dfCast, fetHec < 0.01 & fetEle < 0.01 & fetPar < 0.01)

# ##
# rfam = unique(df$Rfam)[grep("LINE/L2", unique(df$Rfam))[1]]
# ggplot(subset(df, Rfam == rfam)) +
#   geom_bar(aes(x=Rfam, y=Prop*100, group=set, fill=Region, col=Region), stat="identity", position="dodge", alpha=0.6) +
#   facet_wrap(~Species) +
#   ggtitle(rfam)

14/5187*100 # ref
85/7227*100 # spp


###### Enrichment in Repeat Region vs. Background #########
#################### Fisher exact test ####################
#// this tests if the number of repeats of a given family increases in the repeat region as compared to the rest of the genome

# dfCast = dcast(df, Rfam ~ set, value.var = "Count")
# dfCast[is.na(dfCast)] = 0

dfCast$fet2Mel = 0
dfCast$fet2Hec = 0
dfCast$fet2Ele = 0
dfCast$fet2Par = 0

for (rw in 1:nrow(dfCast)) {
  # rfam
  print(dfCast$Rfam[rw])
  # counts specific rfam
  melRRspe = dfCast$hmelv25.RR[rw]; melWGspe = dfCast$hmelv25.WG[rw]
  hecRRspe = dfCast$hhec.RR[rw]; hecWGspe = dfCast$hhec.WG[rw]
  eleRRspe = dfCast$hele.RR[rw]; eleWGspe = dfCast$hele.WG[rw]
  parRRspe = dfCast$hpar.RR[rw]; parWGspe = dfCast$hpar.WG[rw]
  # counts specific rfam
  melRRall = sum(dfCast$hmelv25.RR); melWGall = sum(dfCast$hmelv25.WG)
  hecRRall = sum(dfCast$hhec.RR); hecWGall = sum(dfCast$hhec.WG)
  eleRRall = sum(dfCast$hele.RR); eleWGall = sum(dfCast$hele.WG)
  parRRall = sum(dfCast$hpar.RR); parWGall = sum(dfCast$hpar.WG)
  # matrices for FET
  mat_mel = matrix(c(melWGall,melRRall,melWGspe,melRRspe), ncol=2, dimnames = list(c("back", "target"), c("wg", "rr")))
  mat_hec = matrix(c(hecWGall,hecRRall,hecWGspe,hecRRspe), ncol=2, dimnames = list(c("back", "target"), c("wg", "rr")))
  mat_ele = matrix(c(eleWGall,eleRRall,eleWGspe,eleRRspe), ncol=2, dimnames = list(c("back", "target"), c("wg", "rr")))
  mat_par = matrix(c(parWGall,parRRall,parWGspe,parRRspe), ncol=2, dimnames = list(c("back", "target"), c("wg", "rr")))
  # FET
  fet_mel = fisher.test(mat_hec, alternative = "greater")
  fet_hec = fisher.test(mat_hec, alternative = "greater")
  fet_ele = fisher.test(mat_ele, alternative = "greater")
  fet_par = fisher.test(mat_par, alternative = "greater")
  # add info to table
  dfCast$fet2Mel[rw] = fet_mel$p.value 
  dfCast$fet2Hec[rw] = fet_hec$p.value 
  dfCast$fet2Ele[rw] = fet_ele$p.value 
  dfCast$fet2Par[rw] = fet_par$p.value 
}

nrow(subset(dfCast, fetMel >= 0.01 & fetHec < 0.01 & fetEle < 0.01 & fetPar < 0.01))
subset(dfCast, fetMel >= 0.01 & fetHec < 0.01 & fetEle < 0.01 & fetPar < 0.01)

subset(dfCast, fetHec < 0.01 & fetEle < 0.01 & fetPar < 0.01)

subset(dfCast, 
       fetHec < 0.01 & fetEle < 0.01 & fetPar < 0.01 &
       fet2Mel >= 0.01 & fet2Hec < 0.01 & fet2Ele < 0.01 & fet2Par < 0.01 
)
