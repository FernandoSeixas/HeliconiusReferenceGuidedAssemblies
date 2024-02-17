
##### read data //////////
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/quibl/"
set = "10kb_50kb"
minSites = 500
dat=paste0("chromosomes.phyml_bionj.",set,".data.tsv")
tre=paste0("chromosomes.phyml_bionj.",set,".trees")
data = read.table(paste0(dir,dat))
trees = read.table(paste0(dir,tre))

names(data) = c("scaffold","sta","end","midSites", "sites", "-lnL")

# plot(density(data$sites/10000))
# nrow(subset(data, sites >= 500))/nrow(data)*100

data$index = seq(1,nrow(data),1)
trees$index = seq(1,nrow(trees),1)
validTreeIndex = subset(data, sites >= minSites)$index

dataFilter = subset(data, index %in% validTreeIndex)
treesFilter = subset(trees, index %in% validTreeIndex)

write.table(file=paste0(dir,"chromosomes.phyml_bionj.",set,".minSites",minSites,".trees"), treesFilter[,1], quote=F, row.names = F, col.names = F)
write.table(file=paste0(dir,"chromosomes.phyml_bionj.",set,".minSites",minSites,".data.tsv"), dataFilter[,1:6], quote=F, row.names = F, col.names = F)
