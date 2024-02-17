## load modules
require(stringr)
require(tidyr)
require(ggplot2)
require(reshape2)


##### read data ==================================================
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/TEcounts/"

lfiles = list.files(path = dir, pattern = ".Rfam.")
lfiles = lfiles[!grepl(pattern = "RR02", lfiles)]
# lfiles = lfiles[!grepl(pattern = "RR08", lfiles)]
lfiles = lfiles[!grepl(pattern = "RR09", lfiles)]

df = data.frame() 
for (fl in lfiles) {
  inp = read.table(paste0(dir,fl))
  # inp$prop = inp$V2/sum(inp$V2)
  inp$set = fl
  df = rbind(df, inp)
}
names(df) = c("Rfam","Count","set")
head(df)

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

dfCast$hmelv25.WG = dfCast$hmelv25.WG - dfCast$hmelv25.RR08
dfCast$hhec.WG = dfCast$hhec.WG - dfCast$hhec.RR08
dfCast$hele.WG = dfCast$hele.WG - dfCast$hele.RR08
dfCast$hpar.WG = dfCast$hpar.WG - dfCast$hpar.RR08

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


#################### Enrichment - Fisher Exact Test #########################
dfCast = dcast(df, Rfam ~ set, value.var = "Count")
dfCast[is.na(dfCast)] = 0

# counts all classess
r2 = sum(dfCast$hmelv25.RR); R2 = sum(dfCast$hmelv25.WG);
h2 = sum(dfCast$hhec.RR); H2 = sum(dfCast$hhec.WG)
e2 = sum(dfCast$hele.RR); E2 = sum(dfCast$hele.WG)
p2 = sum(dfCast$hpar.RR); P2 = sum(dfCast$hpar.WG)

# Overall Proportion of repeats in the RR/WG in the target species (hec/ele/par) versus the reference genome (hmelv25)
(h2/H2)/(r2/R2)
(e2/E2)/(r2/R2)
(p2/P2)/(r2/R2)

mat_hec_all = matrix(c(R2,r2,H2,h2), ncol=2, dimnames = list(c("wg", "rr"),c("ref", "hec")))
mat_ele_all = matrix(c(R2,r2,E2,e2), ncol=2, dimnames = list(c("wg", "rr"),c("ref", "ele")))
mat_par_all = matrix(c(R2,r2,P2,p2), ncol=2, dimnames = list(c("wg", "rr"),c("ref", "par")))
fisher.test(mat_hec_all, alternative = "greater")$p.value
fisher.test(mat_ele_all, alternative = "greater")$p.value
fisher.test(mat_par_all, alternative = "greater")$p.value

##### Enrichment of particular Rfam (Repeat families)

# Fisher exact test for enrichment of Rfam/all in RR in relation to WG within species
dfCast$fet1_mel = 0
dfCast$fet1_hec = 0
dfCast$fet1_ele = 0
dfCast$fet1_par = 0
# Fisher exact test for enrichment of Rfam/all in RR in relation to another species
# Q: Is a particular family of repeats more enriched in the RR over the WG in the target species than in the reference species? 
dfCast$fet2_hec = 0
dfCast$fet2_ele = 0
dfCast$fet2_par = 0

rw=21
for (rw in 1:nrow(dfCast)) {
  # rfam
  print(dfCast$Rfam[rw])
  # counts per class
  r1 = dfCast$hmelv25.RR[rw]; R1 = dfCast$hmelv25.WG[rw];
  h1 = dfCast$hhec.RR[rw]; H1 = dfCast$hhec.WG[rw]
  e1 = dfCast$hele.RR[rw]; E1 = dfCast$hele.WG[rw]
  p1 = dfCast$hpar.RR[rw]; P1 = dfCast$hpar.WG[rw]
  # mat 1
  mat_mel1 = matrix(c(R2,R1,r2,r1), ncol=2, dimnames = list(c("all", "rfam"),c("wg", "rr")))
  mat_hec1 = matrix(c(H2,H1,h2,h1), ncol=2, dimnames = list(c("all", "rfam"),c("wg", "rr")))
  mat_ele1 = matrix(c(E2,E1,e2,e1), ncol=2, dimnames = list(c("all", "rfam"),c("wg", "rr")))
  mat_par1 = matrix(c(P2,P1,p2,p1), ncol=2, dimnames = list(c("all", "rfam"),c("wg", "rr")))
  # fet 1
  dfCast$fet1_mel[rw] = fisher.test(mat_mel1, alternative = "greater")$p.value
  dfCast$fet1_hec[rw] = fisher.test(mat_hec1, alternative = "greater")$p.value
  dfCast$fet1_ele[rw] = fisher.test(mat_ele1, alternative = "greater")$p.value
  dfCast$fet1_par[rw] = fisher.test(mat_par1, alternative = "greater")$p.value
  # mat 2
  mat_hec2 = matrix(c(r2,r1,h2,h1), ncol=2, dimnames = list(c("wg", "rr"),c("ref", "hec")))
  mat_ele2 = matrix(c(r2,r1,e2,e1), ncol=2, dimnames = list(c("wg", "rr"),c("ref", "ele")))
  mat_par2 = matrix(c(r2,r1,p2,p1), ncol=2, dimnames = list(c("wg", "rr"),c("ref", "par")))
  # fet 2
  dfCast$fet2_hec[rw] = fisher.test(mat_hec2, alternative = "greater")$p.value
  dfCast$fet2_ele[rw] = fisher.test(mat_ele2, alternative = "greater")$p.value
  dfCast$fet2_par[rw] = fisher.test(mat_par2, alternative = "greater")$p.value
  
}


subset(dfCast, fet1_hec < 0.01 & fet1_ele < 0.01 & fet1_par < 0.01)
subset(dfCast, fet2_hec < 0.01 & fet2_ele < 0.01 & fet2_par < 0.01)

