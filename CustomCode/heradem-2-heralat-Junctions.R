#### load modules =====================================
require(stringr)
require(ggplot2)
require(ggrepel)
require(foreach)
require(doParallel)

#### variables ========================================
filenames = "/n/scratchlfs02/mallet_lab/fseixas/1.pseudo_references/heradem-2-heralat/heradem-2-heralat.mm2"
minAlign = 1000
minMQ = 60
maxGap = 20000
maxDelta = 2.5
minPropMap = 0.5


#### read data and analyze ============================
## read data and apply basic filters [divergence and minimum alignment length] -----
inp=read.table(filenames, fill=T)
names(inp)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
# filter alignments shorter than "minAlign"
inp = subset(inp, Nb.Bases >= minAlign)
# filter alignments too divergent
inp$dv = as.numeric(str_replace(inp$dv, "dv:f:", ""))
inp = subset(inp, dv < 0.10)
# filter low MQ alignments
inp = subset(inp, MQ >= minMQ)

#
inp$distSta = inp$R.start 
inp$distEnd = inp$R.length - inp$R.end 
inp$MinDistTip = 0
inp$Qmap = inp$Q.end-inp$Q.start
# split chromosome and scaffold numbers
inp$chrom = substr(inp$Q.scaffold, 7, 8); inp$chrom = str_remove(inp$chrom, "Hmel2")
inp$scaNB = substr(inp$Q.scaffold, 9, 10)
#
for (i in seq(1,9,1)) {
  oldchr = paste0("_chr", i, "_")
  newchr = paste0("_chr0", i, "_")
  inp$R.scaffold = str_replace(inp$R.scaffold, oldchr, newchr)
}
inp$Qchr = substr(inp$R.scaffold, 8, 9); 
inp$Ssca = substr(inp$R.scaffold, 11, 12); 


###
droplevels(subset(inp, Q.scaffold == "Herato1103"))

aaa = droplevels(subset(inp, R.scaffold == "Hel_chr03_5" & chrom == "03"))
aaa = aaa[order(aaa$Q.start),1:12]
table(aaa$Q.scaffold)
ch = "Herato0309"
sum(subset(aaa, Q.scaffold == ch & strand == "+")$Nb.Bases)/unique(subset(aaa, Q.scaffold == ch)$Q.length)*100
sum(subset(aaa, Q.scaffold == ch & strand == "-")$Nb.Bases)/unique(subset(aaa, Q.scaffold == ch)$Q.length)*100
# subset(aaa, Q.scaffold == ch)

# aaa
# tail(aaa[order(aaa$Q.start),1:12], 100)

aaa = droplevels(subset(inp, Q.scaffold == "Herato0310"))
aaa = aaa[order(aaa$Q.start),1:12]
bbb=aggregate(aaa$Q.end-aaa$Q.start, list(aaa$R.scaffold), sum)
bbb$x = bbb$x/unique(aaa$Q.length)*100
bbb[bbb$x > 1,]
