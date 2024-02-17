#### load modules =====================================
require(stringr)
require(ggplot2)
require(ggrepel)
require(foreach)
require(doParallel)

#### variables ========================================
filenames = "/n/scratchlfs02/mallet_lab/fseixas/1.pseudo_references/heradem-2-heralat/heralat-2-heradem.mm2"
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
inp$chrom = substr(inp$R.scaffold, 7, 8); inp$chrom = str_remove(inp$chrom, "Hmel2")
inp$scaNB = substr(inp$R.scaffold, 9, 10)
for (i in seq(1,9,1)) {
  oldchr = paste0("_chr", i, "_")
  newchr = paste0("_chr0", i, "_")
  inp$Q.scaffold = str_replace(inp$Q.scaffold, oldchr, newchr)
}
inp$Qchr = substr(inp$Q.scaffold, 8, 9); 
inp$Ssca = substr(inp$Q.scaffold, 11, 12); 
head(inp)

# retain big chroms only
bigchrom = c(paste0("0", seq(1,9,1)), seq(10,21,1))
inp = subset(inp, chrom %in% bigchrom)
head(inp)

## retain putative linking contigs
PutativeLinks = data.frame()
for (contig in unique(inp$Q.scaffold)) {
  s1 = subset(inp, Q.scaffold == contig)
  s1 = s1[order(s1$Q.start),]
  # only linking same chromosome
  QChrom = unique(s1$Qchr) 
  s2 = subset(s1, chrom == QChrom)
  # check if maps to multiple scaffolds within one chromosome
  NMChrom = length(unique(s2$chrom))
  NMScaff = length(unique(s2$R.scaffold))
  # only accept scaffolds mapping to one chromosome but different scaffolds
  if (NMChrom == 1 & NMScaff > 1) {}
    for (rw in 1:nrow(s2)) { 
      s2$MinDistTip[rw] = min(s2$distSta[rw], s2$distEnd[rw]) 
      s3 = subset(s2, MinDistTip < 50000)
    }
    if (nrow(s3) > 1 & length(unique(s3$R.scaffold)) > 1) { PutativeLinks = rbind(PutativeLinks, s3) }
  }
}

head(PutativeLinks)
nrow(PutativeLinks)

PutativeLinksSimp = PutativeLinks[,c(1,2,3,4,21,6,22,23,7,8,9,20,5)]
head(PutativeLinksSimp)

subset(PutativeLinksSimp, chrom == "02" & Q.scaffold == "Hel_chr02_1")

aaa = droplevels(subset(inp, Q.scaffold == "Hel_chr17_6" & chrom == "17"))
aaa = aaa[order(aaa$Q.start),1:12]
table(aaa$R.scaffold)
ch = "Herato1718"
# subset(aaa, R.scaffold == ch)
sum(subset(aaa, R.scaffold == ch & strand == "+")$Nb.Bases)/unique(subset(aaa, R.scaffold == ch)$R.length)*100
sum(subset(aaa, R.scaffold == ch & strand == "-")$Nb.Bases)/unique(subset(aaa, R.scaffold == ch)$R.length)*100

# aaa
# tail(aaa[order(aaa$Q.start),1:12], 100)

aaa = droplevels(subset(inp, R.scaffold == "Herato1802"))
aaa = aaa[order(aaa$Q.start),1:12]
bbb=aggregate(aaa$Nb.Bases, list(aaa$Q.scaffold), sum)
bbb$x = bbb$x/unique(aaa$R.length)*100
bbb[bbb$x > 1,]
