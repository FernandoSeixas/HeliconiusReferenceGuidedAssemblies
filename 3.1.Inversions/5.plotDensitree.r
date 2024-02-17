# load libraries
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("stringr")
library("phangorn")



#################### variables ####################
ref = "hmelv25" 
cla = paste0(ref,"_clade_ext")

dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/3.iqtree/"
lfiles = list.files(dir, pattern = "phyml_bionj.w5000.trees.gz")
lfiles = lfiles[grep(ref, lfiles)]
lfiles = lfiles[grep(cla, lfiles)]

lfiles
fil=lfiles[1]
fil


#################### read data ####################
# read tree
nwk=paste0(dir,fil)
tree <- read.tree(nwk)
length(tree)

# RF distance
# n=nrow(tree[[1]]$edge); 2*(n-3)
mean(as.vector(RF.dist(tree)))
# mean(as.vector(KF.dist(tree)))

# update tip labels
for (tr in 1:length(tree)) {
  # update [normalize] names
  tree[[tr]]$tip.label = ifelse(tree[[tr]]$tip.label == "HEL_ERAtog_HEL_ERAtog_A", "HEL_ERA_HEL_ERA_A", tree[[tr]]$tip.label)
  tree[[tr]]$tip.label = ifelse(tree[[tr]]$tip.label == "HEL_ERAtog_HEL_ERAtog_B", "HEL_ERA_HEL_ERA_B", tree[[tr]]$tip.label)
  tree[[tr]]$tip.label = ifelse(tree[[tr]]$tip.label == "HEL_HIM_fat_HEL_HIM_fat_A", "HEL_HIM_HEL_HIM_A", tree[[tr]]$tip.label)
  tree[[tr]]$tip.label = ifelse(tree[[tr]]$tip.label == "HEL_HIM_fat_HEL_HIM_fat_B", "HEL_HIM_HEL_HIM_B", tree[[tr]]$tip.label)
  tree[[tr]]$tip.label = substr(tree[[tr]]$tip.label, 9, 17)
}

# normalize branch lengths
for (tr in 1:length(tree)) { tree[[tr]] = compute.brlen(tree[[tr]]) }

# define baseTree
if (cla == "hmelv25_clade_res") {baseTree=read.tree(text="(((((HEL_CYD_A,HEL_CYD_B),(HEL_TIM_A,HEL_TIM_B)),(HEL_MEL_A,HEL_MEL_B)),((HEL_HEC_A,HEL_HEC_B),((HEL_ELE_A,HEL_ELE_B),(HEL_PAR_A,HEL_PAR_B)))),((HEL_NUM_A,HEL_NUM_B),(HEL_BES_A,HEL_BES_B)));") }
if (cla == "heradem_clade_res") {baseTree=read.tree(text="(((((HEL_HIM_A,HEL_HIM_B),(HEL_ERA_A,HEL_ERA_B)),(HEL_SIA_A,HEL_SIA_B)),(HEL_TEL_A,HEL_TEL_B)),((HEL_DEM_A,HEL_DEM_B),(HEL_SAR_A,HEL_SAR_B)));") }
if (cla == "hmelv25_clade_ext") {baseTree=read.tree(text="(((((HEL_CYD_A,HEL_CYD_B),(HEL_TIM_A,HEL_TIM_B)),(HEL_MEL_A,HEL_MEL_B)),(((HEL_BES_A,HEL_BES_B),(HEL_NUM_A,HEL_NUM_B)),((HEL_HEC_A,HEL_HEC_B),((HEL_ELE_A,HEL_ELE_B),(HEL_PAR_A,HEL_PAR_B))))),(HEL_BUR_A,HEL_BUR_B));") }
if (cla == "heradem_clade_ext") {baseTree=read.tree(text="((((((HEL_HIM_A,HEL_HIM_B),(HEL_ERA_A,HEL_ERA_B)),(HEL_SIA_A,HEL_SIA_B)),(HEL_TEL_A,HEL_TEL_B)),((HEL_DEM_A,HEL_DEM_B),(HEL_SAR_A,HEL_SAR_B))),(HEL_BUR_A,HEL_BUR_B));") }

# densitry
fil
densiTree(tree, consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)
length(tree)
densiTree(tree[40], consensus = baseTree, type = "cladogram", col="blue", nodes="intermediate", use.edge.length=TRUE)

comparePhylo(tree[1], tree[2], plot=TRUE, use.edge.length = TRUE)

require("Quartet")
QuartetDistance(tree[1], tree[2])
