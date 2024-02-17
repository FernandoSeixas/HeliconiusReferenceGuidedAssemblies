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



#################### Enrichment #########################
dfCast = dcast(df, Rfam ~ set, value.var = "Count")
dfCast[is.na(dfCast)] = 0

# Proportions
dfCast$Prop1_hec = 0
dfCast$Prop1_ele = 0
dfCast$Prop1_par = 0
dfCast$Prop2_hec = 0
dfCast$Prop2_ele = 0
dfCast$Prop2_par = 0

# counts all classess
r2 = sum(dfCast$hmelv25.RR); R2 = sum(dfCast$hmelv25.WG);
h2 = sum(dfCast$hhec.RR); H2 = sum(dfCast$hhec.WG)
e2 = sum(dfCast$hele.RR); E2 = sum(dfCast$hele.WG)
p2 = sum(dfCast$hpar.RR); P2 = sum(dfCast$hpar.WG)

for (rw in 1:nrow(dfCast)) {
  # rfam
  print(dfCast$Rfam[rw])
  # counts per class
  r1 = dfCast$hmelv25.RR[rw]; R1 = dfCast$hmelv25.WG[rw];
  h1 = dfCast$hhec.RR[rw]; H1 = dfCast$hhec.WG[rw]
  e1 = dfCast$hele.RR[rw]; E1 = dfCast$hele.WG[rw]
  p1 = dfCast$hpar.RR[rw]; P1 = dfCast$hpar.WG[rw]
  # Enrichment of specific Rfam (genomic background as null hypothesis)
  dfCast$Prop1_hec[rw] = (h1/h2)/(H1/H2); h1+h2+H1+H2
  dfCast$Prop1_ele[rw] = (e1/e2)/(E1/E2); e1+e2+E1+E2
  dfCast$Prop1_par[rw] = (p1/p2)/(P1/P2); p1+p2+P1+P2
  # Enrichment of specific Rfam (hmelv25 as null hypothesis)
  dfCast$Prop2_hec[rw] = (h1/h2)/(r1/r2); h1+h2+r1+r2 
  dfCast$Prop2_ele[rw] = (e1/e2)/(r1/r2); e1+e2+r1+r2 
  dfCast$Prop2_par[rw] = (p1/p2)/(r1/r2); p1+p2+r1+r2 
}

# remove Na and Inf
dfCast[is.na(dfCast)] = 0
for (cl in 10:ncol(dfCast)) {
  for (rw in 1:nrow(dfCast)) {
    dfCast[rw,cl] = ifelse(is.infinite(dfCast[rw,cl]), 0, dfCast[rw,cl])
  }
}


## Outliers
subset(dfCast,
       Prop1_hec > 1.05 & Prop1_ele > 1.05 & Prop1_par > 1.05 &
       Prop2_hec > 1.05 & Prop2_ele > 1.05 & Prop2_par > 1.05
)

subset(dfCast, Prop2_hec > 1 & Prop2_ele > 1 & Prop2_par > 1 )

# subset(dfCast, Prop1_hec > 1 & Prop1_ele > 1 & Prop1_par > 1 )
# 
# 
# ##
# rfam = unique(df$Rfam)[grep("DNA/PiggyBac", unique(df$Rfam))[1]]
# ggplot(subset(df, Rfam == rfam)) +
#   geom_bar(aes(x=Rfam, y=Prop*100, group=set, fill=Region, col=Region), stat="identity", position="dodge", alpha=0.6) +
#   facet_wrap(~Species) +
#   ggtitle(rfam)

