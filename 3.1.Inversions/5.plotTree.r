# load libraries
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("stringr")
library("phangorn")
require("gridExtra")
require("svglite")


#################### variables ####################
ref = "hmelv25" 
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/CladeInversions2Tree/"
# lfiles = list.files(dir, pattern = "outgroupNer.phy.treefile")
lfiles = list.files(dir, pattern = "outgroup.phy.treefile")
lfiles = lfiles[grep(ref, lfiles)]
lfiles


for (i in 1:length(lfiles)) {
  fil=lfiles[i]
  print(fil)
  
  #################### read data ####################
  # read tree
  nwk=paste0(dir,fil)
  tree <- read.tree(nwk)
  tree$tip.label = substr(tree$tip.label, 1, 7)
  
  # remove outgroup (eueides) for visualization
  rootTaxon = "EUE_TAL"
  noOut = drop.tip(tree, rootTaxon)
  
  # transform tip names
  noOut$tip.label = tolower(noOut$tip.label)
  noOut$tip.label = str_replace(noOut$tip.label, "hel_", "h")
  noOut$tip.label = str_replace(noOut$tip.label, "lap_", "h")

  # cls <- list(c1=c("hmel", "hnum"),
  #             c2=c("hbur", "hdor"),
  #             c3=c("hera", "hsia")
  #             )
  # noOut = groupOTU(noOut, cls)
  # Clades = data.frame(
  #   Species = c("hmel","hnum","hbur","hdor","hera","hsia"),
  #   Clade = c("mel","mel","mid","mid","era","era")
  # )
  # 
  
  #################### plot phylogenies ####################

  # tittle
  tittles = LETTERS[seq(1,10,1)]
  if (ref == "hmelv25") { tittle = tittles[i] }
  if (ref == "heradem") { tittle = tittles[i+5] }

  # plot  
  p <- ggtree(noOut, ladderize=TRUE, alpha=0.5) +
    geom_text2(aes(label=label, subset = !isTip & as.numeric(label) > 50), hjust=1.2, vjust=-1.0) +
    geom_tippoint(aes(color=label), size=2, alpha=0.8) +
    geom_tiplab(aes(color=label), size=4, hjust=-0.1) +
    theme(legend.position = "none") +
    # scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3","#7570b3","#d95f02","#1b9e77")) +
    geom_treescale(x = 0.001) +
    ggtitle(tittle) + xlim_tree(0.06) + 
    theme(
      title = element_text(face="bold"), plot.margin = unit(c(0,0,0,0), "lines")
    )
  p
  if (ref == "hmelv25") { assign(paste0("p", (i)),  p) }
  if (ref == "heradem") { assign(paste0("p", (i+5)),  p) }
}

gFinal = grid.arrange(p1, p2, p3, p4, p5,
                      p6, p7, p8, p9, p10,
                      nrow=2)

ggsave(filename = paste0(dir,"InversionPhylogenies.Outgroup.png"), gFinal, width = 14, height = 5, dpi=300)



