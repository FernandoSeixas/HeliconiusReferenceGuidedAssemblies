## read data ==================================================
chr = "hmel21"
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/"
fil = paste0(chr,".target.exons.txt")
inp = read.table(paste0(dir,fil))
names(inp) = c("gene","exon","scaffold","sta","end","strand")
inp$len = inp$end-inp$sta+1

## ==================================================
df = data.frame()
for (gn in unique(inp$gene)) {
  s1 = subset(inp, gene == gn)
  s1$newlen = 0
  s1$frame = 0
  s1$newlen[1] = s1$len[1]
  s1$frame[1] = 1
  s1$excess = 0
  if (nrow(s1) > 1) {
    for (ex in 2:nrow(s1)) {
      ll = s1$newlen[(ex-1)]
      val = round(((ll/3)-floor(ll/3))*3)
      s1$excess[(ex-1)] = val
      if (val == 1) {s1$frame[ex] = 3; s1$newlen[ex] = s1$len[ex] - 2 }
      if (val == 2) {s1$frame[ex] = 2; s1$newlen[ex] = s1$len[ex] - 1 }
      if (val == 0) {s1$frame[ex] = 1; s1$newlen[ex] = s1$len[ex] - 0 }
    }
  }
  # safe check in case I didn't extract all exons from one gene (thus might not be calculating the excess of bases for that)
  if (nrow(s1) == nrow(s1)) {
    ll = s1$newlen[nrow(s1)]
    val=round(((ll/3)-floor(ll/3))*3)
    s1$excess[nrow(s1)] = val
  }
  # add to df
  df = rbind(df, s1)
}
df$frsta = df$frame
df$frend = df$newlen

write.table(file=paste0(dir,chr,".target.exonsFrame.txt"), df, row.names=F, quote=F)


