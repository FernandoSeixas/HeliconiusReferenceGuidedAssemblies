# load libraries
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("stringr")
library("phangorn")

## read tree
# chr="Herato1301";	bk1=22307500; bk2=22668029
# chr="Herato2101"; bk1=10443345; bk2=10815508

chr="Hmel213001o"; bk1=16917575; bk2=17154514
# chr="Hmel221001o"; bk1=8021528; bk2=8342578

# window="w20s20"
window="w40s40"

dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/phyml/"
lfiles = list.files(dir, pattern = paste0(window,".trees.gz"))
lfiles = lfiles[grep(chr,lfiles)]
tree <- read.tree(paste0(dir,lfiles))
length(tree)

## filter trees
filter = read.table(paste0(dir, chr,".phyml_bionj.",window,".data.tsv"), header=T)
pass = as.numeric(rownames(filter[(filter$lnL <= -2000),]))
# pass = as.numeric(rownames(filter[(filter$sites >= 1000),]))
# tree = tree[pass]
length(tree)

# normalize branch lengths
for (tr in 1:length(tree)) { tree[[tr]] = compute.brlen(tree[[tr]]) }

# update species names
for (tr in 1:length(tree)) { tree[[tr]]$tip.label = substr(tree[[tr]]$tip.label, 1, 7) }

# remove outgroup (eueides) for visualization
rootTaxon = "EUE_TAL"
noOut = tree
for (tr in 1:length(noOut)) { noOut[[tr]] = drop.tip(noOut[[tr]], rootTaxon) }

# base tree
baseTree=read.tree(text="((HEL_ERA,HEL_SIA),(LAP_DOR,HEL_BUR),(HEL_MEL,HEL_NUM));")

# densitree
densiTree(noOut, consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)

densiTree(noOut[1], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[2], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[3], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[4], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)

densiTree(noOut[12], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[13], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[14], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[15], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)

densiTree(noOut[c(22)], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)

length(noOut)

text(x=0, y=7, paste0(filter[i,1],":",filter[i,2],"-",filter[i,3]), pos = 4)

densiTree(noOut[4:6], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[7:9], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)

densiTree(noOut[4], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)





densiTree(noOut[6], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)


densiTree(noOut[22], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
densiTree(noOut[23], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)








