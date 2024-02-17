# load libraries
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("stringr")
library("phangorn")



#################### variables ####################
win = "5kb_50kb"
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/quibl/"
# fil = "chromosomes.phyml_bionj.w05b50.subset2.trees"
fil = paste0("chromosomes.phyml_bionj.",win,".trees")
dat = paste0("chromosomes.phyml_bionj.",win,".data.tsv")

nwk = paste0(dir,fil)
data = paste0(dir,dat)


#################### read data ####################
# read tree
tree <- read.tree(nwk)
class(tree)

# data.tsv
dataphylo = read.table(data)
names(dataphylo) = c("scaffold","start","end","mid","sites","lnL")
dataphylo$index = seq(1:nrow(dataphylo))
# dataFilter = subset(dataphylo, lnL < -1000)
dataFilter = subset(dataphylo, sites >= 500)
nrow(dataFilter)/nrow(dataphylo)*100
# plot(density(dataphylo$sites))
# plot(density(dataphylo$lnL))
# abline(v = -500)
plot(x=dataphylo$sites, y=dataphylo$lnL)

# filter trees
treeFilter = tree[dataFilter$index]

# root trees
treeFilter = root.multiPhylo(phy = treeFilter, "eueA")
# treeFilter = root.multiPhylo(phy = treeFilter, "eueB")

# normalize branch lengths
for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = compute.brlen(treeFilter[[tr]]) }

# remove one outgroup 
# for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = drop.tip(phy = treeFilter[[tr]], tip = "eueA") }
for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = drop.tip(phy = treeFilter[[tr]], tip = "eueB") }
# for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = drop.tip(phy = treeFilter[[tr]], tip = "burA") }
# for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = drop.tip(phy = treeFilter[[tr]], tip = "burB") }
# for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = drop.tip(phy = treeFilter[[tr]], tip = "dorA") }
# for (tr in 1:length(treeFilter)) { treeFilter[[tr]] = drop.tip(phy = treeFilter[[tr]], tip = "dorB") }

# define baseTree
baseTree = read.tree(text=("(((((melB,melA),((burB,burA),(dorB,dorA))),(eraA,eraB))),eueA);"))
# baseTree = read.tree(text=("(((((melB,melA),(dorB,dorA)),(eraA,eraB))),eueA);"))
# baseTree = read.tree(text=("(((((melB,melA),(dorB,dorA)),(eraA,eraB))),eueB);"))
# baseTree = read.tree(text=("(((melB,melA),(eraA,eraB)),eueB);"))
baseTree = collapse.singles(baseTree)

# plot densitry
n = 1/length(treeFilter)*10
n = 0.002

densiTree(treeFilter[sample(1:length(treeFilter), 500)], 
          consensus = baseTree, 
          type = "cladogram", 
          col="blue", 
          nodes="intermediate", 
          alpha = n,
          jitter = list(amount=.1, random=FALSE),
          # scaleX = TRUE,
          use.edge.length=TRUE
)


# RF distance
# n=nrow(tree[[1]]$edge); 2*(n-3)
mean(as.vector(RF.dist(tree)))
# mean(as.vector(KF.dist(tree)))


