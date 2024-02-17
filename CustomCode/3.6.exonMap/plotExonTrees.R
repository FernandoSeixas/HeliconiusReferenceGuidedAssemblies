# load libraries
library("ape")
library("ggplot2")
library("ggtree")
library("stringr")
library("Biostrings")
require("ips")
# library("phangorn")

exon = "hmel09.HMEL006651g1.t1.E10"

# read files
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.3.GenomeExpansion/exonMapping/1.Map2CN/5.frmseq/"
lfiles = list.files(dir, pattern = ".fasta")
file = lfiles[grepl(paste0(exon,".frame"), lfiles)]

# read fasta
alignment = read.dna(file = paste0(dir,file), format = "fasta")
alignLeng = length(alignment)/length(labels(alignment))
ParsInfoSites = pis(alignment, what = "absolute")
print(paste0("Sequences: ", length(labels(alignment))))
print(paste0("Alignment Length: ", alignLeng))
print(paste0("Parsimony Informative Sites: ", ParsInfoSites))

## define functions
f <- function(x) {
  # neighbor-joining
  phylo = nj(dist.dna(x))
  # root
  if (length(phylo$tip.label[grepl("hera", phylo$tip.label)]) > 1) {
    phyloRoot = root(phy = phylo, node = MRCA(phylo, phylo$tip.label[grepl("hera", phylo$tip.label)])) 
  }
  if (length(phylo$tip.label[grepl("hera", phylo$tip.label)]) == 1) {
    eraName = phylo$tip.label[grepl("hera", phylo$tip.label)][1]
    phyloRoot = root.phylo(phylo, outgroup = eraName)
  }
  return(phyloRoot)
}  

# neighbor-joining
phyloRoot = f(alignment)

# bootstraps
bstrees = boot.phylo(phyloRoot, alignment, f, quiet = TRUE, B = 100)

hecInds = phyloRoot$tip.label[grepl("hhec", phyloRoot$tip.label)]
eleInds = phyloRoot$tip.label[grepl("hele", phyloRoot$tip.label)]
parInds = phyloRoot$tip.label[grepl("hpar", phyloRoot$tip.label)]
is.monophyletic(phyloRoot, hecInds)
is.monophyletic(phyloRoot, eleInds)
is.monophyletic(phyloRoot, parInds)

# compute branch lengths
# phyloRoot = compute.brlen(phyloRoot)

# rename tips
phyloRoot$tip.label = sapply(str_split(phyloRoot$tip.label, "\\."), "[", 1)
phyloRoot$tip.label = str_remove(phyloRoot$tip.label, "_R_")
phyloRoot$tip.label = substr(phyloRoot$tip.label, 1, 4)


## plot tree =====
ggtree(phyloRoot, ladderize=TRUE, alpha=0.5) +
  geom_text2(aes(subset = !isTip), label=bstrees, hjust=1.5, vjust=-1) +
  geom_tippoint(aes(color=label), size=2, alpha=0.8) +
  geom_tiplab(aes(color=label), size=3) +
  theme(legend.position = "none") +
  ggtitle(paste0(file,"\n", "Parsimony Informative Sites: ", ParsInfoSites)) + xlim_tree(0.01) 

