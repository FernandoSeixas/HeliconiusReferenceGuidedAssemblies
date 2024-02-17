##
options(scipen = 999)


## read scaffold lengths 
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/quibl/"
fil = "hmelv25.anchorScaffoldLengths.txt"
scaLens = read.table(paste0(dir,fil))
names(scaLens) = c("scaffold", "len")
scaLens = subset(scaLens, len > 1000000)


## define windows
wsize = 5000
wstep = 50000

for (ch in unique(scaLens$scaffold)) {
  chrDataFrame = data.frame()
  s1 = subset(scaLens, scaffold == ch)
  st = wstep + 1
  en = st + wsize - 1
  mxlen = s1$len
  # df = data.frame(sca = ch, sta = st, end = en)
  # chrDataFrame = rbind(chrDataFrame, df)
  while(en < mxlen) {
    # add to data frame
    df = data.frame(sca = ch, sta = st, end = en)
    chrDataFrame = rbind(chrDataFrame, df)
    # update coordinates
    st = en + wstep + 1
    en = st + wsize - 1
  }
  write.table(file=paste0(dir, "support/", ch, ".phyml.coord.",wsize/1000,"kb_",wstep/1000,"kb.txt"), chrDataFrame, row.names = F, col.names = F, quote=F, sep="\t")
}

