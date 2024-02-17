## load modules
require(ggplot2)
require(ggrepel)
library(gggenes)
require(stringr)
require(gridExtra)

## define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")
tal_clade = c("etal")
refs = c("hmelv25","heradem","heralat")


## subset inversions regions
ref="heradem"

if (ref == "hmelv25") {
  sc="Hmel202001o"; bk1=2147755; bk2=2230606
  # sc="Hmel206001o"; bk1=12129405; bk2=12591271; bk3=12279100; bk4=12280171
  # sc="Hmel213001o"; bk1=16920997; bk2=17154514
  # sc="Hmel221001o"; bk1=8029417; bk2=8342578
}
if (ref == "heradem") {
  sc="Herato0206"; bk1=880283; bk2=1059137
  # sc="Herato0606"; bk1=12189406; bk2=12756548; bk3=12361782; bk4=12367112
  # sc="Herato1301"; bk1=22307500; bk2=22668029
  # sc="Herato2101"; bk1=10443345; bk2=10815508
}

# ref chromosome coordinates
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/"
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)
sb = subset(scaffLengths, scaffold == as.character(sc))
st = sb$c_sta 
en = sb$c_end

print(paste0(sc,":",bk1,"-",bk2))

#################### read original minimap2 mappings ####################
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/"
subdir = "2.minimap2"
filenames <- list.files(paste0(dir, subdir), pattern=paste0("*",ref,".filter.mm2"), full.names=TRUE)
# read data
# if (!exists("minimaps")) {
  minimaps = data.frame()
  for (file in filenames) {
    # define species
    spp = str_replace(file, paste0(dir, subdir,"/"), "")
    spp = str_replace(spp, paste0("-2-",ref,".filter.mm2"), "")
    print(spp)
    # read file
    inp=read.table(paste0(dir, subdir, "/", spp, "-2-",ref,".filter.mm2"))
    names(inp)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
    # add species
    inp$species = spp
    # add to dataframe
    minimaps = rbind(minimaps, inp)
  }
# }


## get scaffolds within specified region
bf = 100000
scaffoldRegions = data.frame()
sub = subset(minimaps, R.scaffold == as.character(sc))
for (rw in 1:nrow(sub)) {
  stSub = sub$R.start[rw]
  enSub = sub$R.end[rw]
  if (stSub > (bk1-bf) & stSub < (bk2+bf) | enSub > (bk1-bf) & enSub < (bk2+bf)) {scaffoldRegions = rbind(scaffoldRegions, sub[rw,])}
}



## adjust start and end depending on strand
for (rw in 1:nrow(scaffoldRegions)) {
  if (scaffoldRegions$strand[rw] == "-") {
    tmp = scaffoldRegions$R.start[rw]
    scaffoldRegions$R.start[rw] = scaffoldRegions$R.end[rw]
    scaffoldRegions$R.end[rw] = tmp
  }
}

# define species order
scaffoldRegions$species = factor(scaffoldRegions$species, levels = c(tal_clade,mel_clade,mid_clade,era_clade,refs))

# define if scaffold maps to 1 (normal) or 2 (inversion) strands
scaffoldRegions$type = 1
scaffoldRegions$span = -1
for (spp in levels(scaffoldRegions$species) ) {
  print(spp)
  # go through each species 
  s1 = subset(scaffoldRegions, species == spp)
  for (contig in unique(s1$Q.scaffold)) {
    # go through each Q.scaffold 
    s2 = subset(s1, Q.scaffold == as.character(contig))
    # must be a split alignment
    if (length(unique(s2$strand)) == 2) {
      s3 = subset(s2, (R.start > bk1 & R.start < bk2) | (R.end > bk1 & R.end < bk2) )
      # check if split aligment is in inversion
      if (nrow(s3) > 0) {
        scaffoldRegions$span = ifelse(scaffoldRegions$species == spp & scaffoldRegions$Q.scaffold == as.character(contig), abs(min(s2$R.start,s2$R.end)-max(s2$R.start,s2$R.end)), scaffoldRegions$span)
        scaffoldRegions$type = ifelse(scaffoldRegions$species == spp & scaffoldRegions$Q.scaffold == as.character(contig), 2, scaffoldRegions$type)
      }
      # print(s2)
    }
  }
}
scaffoldRegions$type = as.character(scaffoldRegions$type)


## get scaffolds overlapping breakpoints
scaffoldBreaks = data.frame()
sub = subset(minimaps, R.scaffold == as.character(sc) & R.start > (bk1-bf) & R.end < (bk2+bf) & MQ >= 10)
sub$breakOverlap = -1
sub$breakType = NA
for (rw in 1:nrow(sub)) {
  stSub = sub$R.start[rw]
  enSub = sub$R.end[rw]
  # overlap breakpoint | straight
  if (bk1 > stSub & bk1 < enSub) { ov = min(abs(bk1-stSub),abs(bk1-enSub)); sub$breakOverlap[rw] = ov; sub$breakType[rw] = "straight"; scaffoldBreaks = rbind(scaffoldBreaks, sub[rw,]) }
  if (bk2 > stSub & bk2 < enSub) { ov = min(abs(bk2-stSub),abs(bk2-enSub)); sub$breakOverlap[rw] = ov; sub$breakType[rw] = "straight"; scaffoldBreaks = rbind(scaffoldBreaks, sub[rw,]) }
  # if (bk3 > stSub & bk3 < enSub) { ov = min(abs(bk3-stSub),abs(bk3-enSub)); sub$breakOverlap[rw] = ov; sub$breakType[rw] = "straight"; scaffoldBreaks = rbind(scaffoldBreaks, sub[rw,]) }
  # if (bk4 > stSub & bk4 < enSub) { ov = min(abs(bk4-stSub),abs(bk4-enSub)); sub$breakOverlap[rw] = ov; sub$breakType[rw] = "straight"; scaffoldBreaks = rbind(scaffoldBreaks, sub[rw,]) }
}

# get rid of very small overlaps [< 1000 bp ]
scaffoldBreaks = subset(scaffoldBreaks, breakOverlap >= 1000)

## adjust scaffoldBreaks start and end coordinates depending on strand
for (rw in 1:nrow(scaffoldBreaks)) {
  if (scaffoldBreaks$strand[rw] == "-") {
    tmp = scaffoldBreaks$R.start[rw]
    scaffoldBreaks$R.start[rw] = scaffoldBreaks$R.end[rw]
    scaffoldBreaks$R.end[rw] = tmp
  }
}

# define species order
scaffoldBreaks$species = factor(scaffoldBreaks$species, levels = c(tal_clade,mel_clade,mid_clade,era_clade,refs))

# filter alignments based on MQ
scaffoldBreaksFilter = subset(scaffoldBreaks, MQ >= 10)
scaffoldRegionsFilter = subset(scaffoldRegions, MQ >= 10)

# plot
if (ref == "hmelv25") {
  a = ggplot() +
    # add breakpoints
    geom_vline(xintercept = bk1/1000000, col="red") +
    geom_vline(xintercept = bk2/1000000, col="red") +
    # geom_vline(xintercept = bk3/1000000, col="red") +
    # geom_vline(xintercept = bk4/1000000, col="red") +
    # plot DISCOVAR scaffolds
    geom_gene_arrow(data=subset(scaffoldRegionsFilter, species %in% c("etal")), aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
    # plot DISCOVAR scaffolds
    geom_gene_arrow(data=subset(scaffoldRegionsFilter, type==2), aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
    scale_fill_manual(values = c("#fdae61","#4e7dff")) +
    scale_color_manual(values = c("#fdae61","#4e7dff")) +
    # scaffolds overlapping breakpoints
    geom_gene_arrow(data=scaffoldBreaksFilter, aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold), fill="black", alpha=1.0) +
    # visual
    xlab("Chromosome Position (Mb)") + ylab("DISCOVAR/w2rap Scaffold") +
    facet_grid(species ~ ., drop = TRUE, scales = "free_y") +
    ggtitle(paste0(ref," | ", sc,":",bk1,"-",bk2)) +
    theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))
    a
    ggsave(filename = paste0(dir,sc,".",bk1,".",bk2,".",ref,".eps"), plot = a, dpi = "print")
}

if (ref == "heradem") {
  b = ggplot() +
    # add breakpoints
    geom_vline(xintercept = bk1/1000000, col="red") +
    geom_vline(xintercept = bk2/1000000, col="red") +
    # geom_vline(xintercept = bk3/1000000, col="red") +
    # geom_vline(xintercept = bk4/1000000, col="red") +
    # plot DISCOVAR scaffolds
    geom_gene_arrow(data=subset(scaffoldRegionsFilter, species == "etal"), aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
    # plot DISCOVAR scaffolds
    geom_gene_arrow(data=subset(scaffoldRegionsFilter, type==2), aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
    scale_fill_manual(values = c("#fdae61","#4e7dff")) +
    scale_color_manual(values = c("#fdae61","#4e7dff")) +
    # scaffolds overlapping breakpoints
    geom_gene_arrow(data=scaffoldBreaksFilter, aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold), fill="black", alpha=1.0) +
    # visual
    xlab("Chromosome Position (Mb)") + ylab("DISCOVAR/w2rap Scaffold") +
    facet_grid(species ~ ., drop = TRUE, scales = "free_y") +
    ggtitle(paste0(ref," | ", sc,":",bk1,"-",bk2)) +
    theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))
    b
    ggsave(filename = paste0(dir,sc,".",bk1,".",bk2,".",ref,".eps"), plot = b, dpi = "print")
}

comb = grid.arrange(b, a, ncol=2)
ggsave(filename = paste0(dir,sc,".inversions.comb.eps"), plot = comb, dpi = "print", scale = 1.5)

# 
# ###
# ggplot() +
#   # add breakpoints
#   geom_vline(xintercept = bk1/1000000, col="red") +
#   geom_vline(xintercept = bk2/1000000, col="red") +
#   # plot DISCOVAR scaffolds
#   geom_gene_arrow(data=subset(scaffoldRegionsFilter, species == "etal"), aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
#   # color
#   scale_fill_manual(values = c("#fdae61","#4e7dff")) +
#   scale_color_manual(values = c("#fdae61","#4e7dff")) +
#   # visual
#   xlab("Chromosome Position (Mb)") + ylab("DISCOVAR/w2rap Scaffold") +
#   facet_grid(species ~ ., drop = TRUE, scales = "free_y") +
#   ggtitle(paste0(ref," | ", sc,":",bk1,"-",bk2)) +
#   theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))
# ###
# subset(scaffoldRegionsFilter, type==2 & species == "hbur" & strand == "+")


unique(subset(sub, species == "hdor")$Q.scaffold)
unique(subset(sub, species == "hbur")$Q.scaffold)

mmm = subset(sub, species == "hdor")
mmm[order(mmm$R.start),]

# mmm = subset(scaffoldBreaksFilter, species == "hnum")
mmm = subset(scaffoldRegionsFilter, type==1 & species == "hnum")

mmm = subset(scaffoldRegionsFilter, type==2 & species == "hdor")
mmm = subset(scaffoldRegionsFilter, type==2 & species == "hbur" & Q.scaffold == "xfSc0000958")
mmm = subset(scaffoldRegionsFilter, type==2 & species == "hbur" & Q.scaffold == "xfSc0000796")
# subset(scaffoldBreaks, species == "hsia")
mmm[order(mmm$R.start),]
# mmm[order(mmm$Q.start),]

sum(mmm$M.Bases)/mmm$Q.length[1]*100

# table(droplevels(mmm$Q.scaffold))





# # plot
# ggplot(scaffoldBreaks) +
#   # plot scaffolds
#   geom_gene_arrow(data=subset(scaffoldRegions, type==2), aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
#   # plot scaffolds
#   geom_gene_arrow(aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold)) +
#   # add breakpoints
#   geom_vline(xintercept = bk1/1000000, alpha=0.8,) +
#   geom_vline(xintercept = bk2/1000000, alpha=0.8,) +
#   # geom_vline(xintercept = bk3/1000000, alpha=0.8) +
#   # geom_vline(xintercept = bk4/1000000, alpha=0.8) +
#   # visual
#   facet_wrap(~species, scales="free_y", ncol=2) +
#   ggtitle(paste0(ref," // ", sc,":",bk1,"-",bk2))





# # plot
# # ggplot(scaffoldRegions, aes(label=species, x=(R.start+R.end)/2000000, y=species, fill=type, colour=type)) +
# ggplot(scaffoldRegions) +
#   # geom_gene_arrow(aes(xmin = R.start/1000000, xmax = R.end/1000000, y = species), alpha=ifelse(scaffoldRegions$type == 2, 1, 0.1)) +
#   geom_gene_arrow(aes(xmin = R.start/1000000, xmax = R.end/1000000, y = species), fill=ifelse(scaffoldRegions$strand == "-", "black", "white")) +
#   geom_vline(xintercept = bk1/1000000, alpha=0.8) +
#   geom_vline(xintercept = bk2/1000000, alpha=0.8)
# 
# ggplot(subset(scaffoldRegions, type==2)) +
#   geom_gene_arrow(aes(xmin = R.start/1000000, xmax = R.end/1000000, y = Q.scaffold, fill=strand, color=strand)) +
#   geom_vline(xintercept = bk1/1000000, alpha=0.8) +
#   geom_vline(xintercept = bk2/1000000, alpha=0.8) +
#   # geom_vline(xintercept = 9986426/1000000, alpha=0.8) +
#   # geom_vline(xintercept = 10185405/1000000, alpha=0.8) +
#   facet_wrap(~species, scales="free_y", ncol=2) +
#   ggtitle(paste0(ref," // ", sc,":",bk1,"-",bk2))