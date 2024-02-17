## load modules
require(ggplot2)
require(gridExtra)

########################## read data  ############################
ref = "hmelv25"
spp = "lapdor"
tre = "tree1"
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/fdM/4.dstats/"
fil = paste0("wgenome.",ref,".outgroup.",spp,".1Mb.",tre,".csv")
inp = read.csv(paste0(dir,fil))
inpFilter = subset(inp, sitesUsed >= 1000)

inpFilter$dummypos = seq(1,nrow(inpFilter),1)
ggplot(inpFilter) + geom_point((aes(x=dummypos, y=D))) + ylim(-1,1)

# ggplot(inpFilter) + geom_density(aes(fdM)) + xlim(-1,1)

ggplot(inpFilter) + geom_density((aes(D))) + xlim(-1,1)


######################### Block Jackknife #########################
# load functions
source("~/software/popgen/genomics_general/jackknife.R")

D.stat = function(dataframe) (sum(dataframe$ABBA) - sum(dataframe$BABA)) / (sum(dataframe$ABBA) + sum(dataframe$BABA))

## info about blocks
blocks = inpFilter[,1:3]
n_blocks = nrow(inpFilter)

##
Dstat.jackknife = c()
for (i in 1:n_blocks) {
  s1 = inpFilter[seq(1:n_blocks)[!(seq(1:n_blocks) %in% i)],]
  Dstat.jackknife = c(Dstat.jackknife, D.stat(s1))
}

Dstat.jackknife = rnorm(n_blocks)

D = D.stat(inpFilter)
D

D_mn = mean(Dstat.jackknife)
D_sd = sd(Dstat.jackknife)
D_err <- D_sd/sqrt(n_blocks)

D_Z <- D / D_err
D_p <- 2*pnorm(-abs(D_Z))
D_p
