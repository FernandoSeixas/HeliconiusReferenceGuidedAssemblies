## load modules =====================================
require(stringr)
require(ggplot2)
require(ggrepel)
require(foreach)
require(doParallel)

## variables ========================================
ref="heralat"
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
subdir = "2.minimap2"
filenames <- list.files(paste0(dir, subdir), pattern=paste0("*",ref,".filter.mm2"), full.names=TRUE)
# filenames <- list.files(paste0(dir, "2.minimap2/vOrigin"), pattern=paste0("*",ref,".filter.mm2"), full.names=TRUE)
# filenames <- list.files(paste0(dir, "2.minimap2/vGlued"), pattern=paste0("*",ref,"Glued.filter.mm2"), full.names=TRUE)
minAlign = 1000            ## minimum length of an alignment to be considered 
minMQ = 20                 ## minimum mapping quality to consider an alignment
maxGap = 50000             ## maximum distance between two alignments of the same Q.scaffold to consider concatenation
maxOverlap = 5000          ## maximum overlap between two mappings on different strands
# maxDelta = 2.5             ## 
minPropMap = 0.2           ## minumim proportion of a query scaffold, so that the scaffold might be considered
minAllowedDistance = 50000 ## minumum allowed distance to reference scaffold boundaries

# define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("heralat", "heradem", "hera","hhimfat","hhim","hsia","htel","hdem","hsar")

## DEFINE FUNCTIONS ================================================
source("~/code/heliconius_seixas/1.pseudo_references/4.1.Inversions/2.inversionFunctions.R")



#################### Detect candidate Inversions ####################
# cl <- makePSOCKcluster(4)
registerDoParallel(4)
stime <- Sys.time()
ParsedFiles = foreach(i=1:length(filenames), .combine=rbind, .packages='stringr') %dopar% { ScaffoldParser(filenames[i]) }
etime <- Sys.time()
print(etime-stime)

# subset scaffolds based on whether they are marked as passing or failing filters 
InvCandidates.Pass1 = subset(ParsedFiles, class == "Pass")
InvCandidates.Fail1 = subset(ParsedFiles, class == "Fail")

# write scaffolds passing and failing filters to different files 
write.table(file=paste0(dir,"InvCandidates.Pass1-2-",ref,".txt"), InvCandidates.Pass1, quote=F, row.names=F, sep="\t")
write.table(file=paste0(dir,"InvCandidates.Fail1-2-",ref,".txt"), InvCandidates.Fail1, quote=F, row.names=F, sep="\t")
nrow(InvCandidates.Pass1)
nrow(InvCandidates.Fail1)


#################### Merge close alignment #########################
InvCandidates.Pass1 = read.table(paste0("/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/InvCandidates.Pass1-2-",ref,".txt"), header=T)

InvCandidates.MergeRaw = data.frame()

for (spp in unique(InvCandidates.Pass1$species)) {
  print(spp)
  # analyze one species at the time
  sub1 = subset(InvCandidates.Pass1, species == spp)
  for (contig in unique(sub1$Q.scaffold)) {
    # analyze one scaffold at the time
    sub2 = subset(sub1, Q.scaffold == contig)
    if (length(unique(sub2$R.scaffold)) > 1) {
      sub3 = aggregate(sub2$M.Bases, list(sub2$R.scaffold), sum)
      sub3Chrom = sub3[which(sub3$x == max(sub3$x)),]$Group.1
      sub2 = subset(sub2, R.scaffold == sub3Chrom)
    }
    # sub2 = ooo[order(ooo$Q.start),] 
    plus =  subset(sub2, strand == "+")
    minus = subset(sub2, strand == "-")
    #
    pMerge = plus
    mMerge = minus
    # merge overlapping alignments
    if (nrow(plus) > 1)  {pMerge = mergeRegions(plus) }
    if (nrow(minus) > 1) {mMerge = mergeRegions(minus) }
    if (!is.null(pMerge) & !is.null(mMerge)) {
      # 2nd, check proportion of scaffold that is mapped
      PropMap = (sum(pMerge$Nb.Bases)+sum(mMerge$Nb.Bases))/unique(sub2$Q.length)
      pMerge$propMap = PropMap
      mMerge$propMap = PropMap
      # if still more than 1 mapping per strand select only the largest alignment
      if (nrow(pMerge) > 1 | nrow(mMerge) > 1) {
        pMerge = pMerge[which(pMerge$Nb.Bases == max(pMerge$Nb.Bases)),]
        mMerge = mMerge[which(mMerge$Nb.Bases == max(mMerge$Nb.Bases)),]
      }
      ## 3rd, check if reciprocal match in Query (overlap)
      Qint = 0
      x1 = pMerge$Q.start; x2 = pMerge$Q.end; lenx = x2-x1+1
      y1 = mMerge$Q.start; y2 = mMerge$Q.end; leny = y2-y1+1
      # test if segments overlap
      if (x1 <= y2 & y1 <= x2) {
        Qint = 1 
        Qoverlap = min(x2, y2)-max(x1, y1)+1 # length of overlap
        # if one segment is NOT completely within the other 
        if (lenx != Qoverlap & leny != Qoverlap & Qoverlap < maxOverlap) { Qint = 0}
      }
      ## 4th, check if reciprocal match in Reference  (overlap)
      Rint = 0
      x1 = pMerge$R.start; x2 = pMerge$R.end; lenx = x2-x1+1
      y1 = mMerge$R.start; y2 = mMerge$R.end; leny = y2-y1+1
      # test if segments overlap
      if (x1 <= y2 & y1 <= x2) {
        Rint = 1 
        Roverlap = min(x2, y2)-max(x1, y1)+1 # length of overlap
        # if one segment is NOT completely within the other 
        if (lenx != Roverlap & leny != Roverlap & Roverlap < maxOverlap) { Rint = 0}
      }
      # add to table if passing all filters
      if (Qint == 0 & Rint == 0) {
        mrgScaff = rbind(pMerge,mMerge)
        InvCandidates.MergeRaw = rbind(InvCandidates.MergeRaw, mrgScaff)
      }
    }
  }
}
# update species names
InvCandidates.MergeRaw$species = str_remove(InvCandidates.MergeRaw$species, paste0("-2-",ref,".filter.mm2"))
# filter based on proportion of the scaffod mapped
ggplot(InvCandidates.MergeRaw, aes(propMap)) + geom_histogram(breaks = seq(0,2,0.01))
InvCandidates.Merge1 = subset(InvCandidates.MergeRaw, propMap >= minPropMap)
# InvCandidates.Merge$pIdent = InvCandidates.Merge$M.Bases/abs(InvCandidates.Merge$Q.end-InvCandidates.Merge$Q.start)
# ggplot(InvCandidates.Merge) + geom_density(aes(pIdent))
# InvCandidates.Merge = subset(InvCandidates.Merge, pIdent >= 0.25)

# filter based on the length mapped to each strand
InvCandidates.Merge2 = data.frame()
for (sp in unique(InvCandidates.Merge1$species)) {
  s1 = subset(InvCandidates.Merge1, species == sp)
  # shortMaps = subset(s1, abs(Q.end-Q.start) < 1000)$Q.scaffold
  shortMaps = subset(s1, Nb.Bases < 1000)$Q.scaffold
  length(shortMaps)
  s2 = s1[!(s1$Q.scaffold %in% shortMaps),]
  InvCandidates.Merge2 = rbind(InvCandidates.Merge2, s2)
}
nrow(InvCandidates.Merge2)/nrow(InvCandidates.Merge1)
# write to file
write.table(file=paste0(dir,"InvCandidates.Merge2-2-",ref,".txt"), InvCandidates.Merge2, quote=F, row.names=F, sep="\t")



#################### Define Breakpoints ####################
# create dataframe
InvCandidates.Breakpoints = data.frame()
# parse
for (spp in unique(InvCandidates.Merge2$species)) {
  sub1 = subset(InvCandidates.Merge2, species == spp)
  for (contig in unique(sub1$Q.scaffold)) {
    sub2 = subset(sub1, Q.scaffold == contig)
    Bsta = 0
    Bend = 0
    if (nrow(sub2) == 2) {
      plus  = subset(sub2, strand == "+")
      minus = subset(sub2, strand == "-")
      PQsta = plus$Q.start
      MQsta = minus$Q.start
      PRsta = plus$R.start
      MRsta = minus$R.start
      if (PQsta < MQsta & PRsta < MRsta) {Bsta = plus$R.end;  Bend = minus$R.end }
      if (PQsta > MQsta & PRsta < MRsta) {Bsta = plus$R.sta;  Bend = minus$R.sta }
      if (PQsta < MQsta & PRsta > MRsta) {Bsta = minus$R.end; Bend = plus$R.end  }
      if (PQsta > MQsta & PRsta > MRsta) {Bsta = minus$R.sta; Bend = plus$R.sta  }
      df = data.frame(species = spp, Q.scaffold = contig, R.scaffold = unique(sub2$R.scaffold), Bsta = Bsta, Bend = Bend )
      InvCandidates.Breakpoints = rbind(InvCandidates.Breakpoints, df)
    }
    # if (nrow(sub2) >= 3) {
    #   Bsta = min(sub2$R.start, sub2$R.end)
    #   Bend = max(sub2$R.start, sub2$R.end)
    #   df = data.frame(species = spp, Q.scaffold = contig, R.scaffold = unique(sub2$R.scaffold), Bsta = Bsta, Bend = Bend )
    #   InvCandidates.Breakpoints2 = rbind(InvCandidates.Breakpoints2, df)
    # }
  }
}
# write to file
write.table(file=paste0(dir,"InvCandidates.Breakpoints-2-",ref,".txt"), InvCandidates.Breakpoints, quote=F, row.names=F, sep="\t")



# #################### Candidate Inversions overlapped by REFERENCE inversions ####################
# InvCandidates.BreakpointsFilter = data.frame()
# for (spp in unique(InvCandidates.Breakpoints$species)) {
#   sub1 = subset(InvCandidates.Breakpoints, species == spp)
#   if (ref == "hmelv25") {Control = subset(InvCandidates.Breakpoints, species == "hmel") }
#   if (ref == "heradem") {Control = subset(InvCandidates.Breakpoints, species == "hera") }
#   if (ref == "heralat") {Control = subset(InvCandidates.Breakpoints, species == "hera") }
#   Control = subset(Control, Bend-Bsta < 3000000)
#   for (contig in unique(sub1$Q.scaffold)) {
#     sub2 = subset(sub1, Q.scaffold == contig)
#     god2 = subset(Control, R.scaffold == as.character(sub2$R.scaffold))
#     fail = 0
#     if (nrow(god2) >= 1) {
#       for (rw in 1:nrow(god2)) {
#         Gsta = god2$Bsta[rw]
#         Gend = god2$Bend[rw]
#         if(sub2$Bsta >= Gsta & sub2$Bsta <= Gend) {fail = 1}
#         if(sub2$Bsta >= Gsta & sub2$Bsta <= Gend) {fail = 1}
#       }
#     }
#     sub2$fail = fail
#     InvCandidates.BreakpointsFilter = rbind(InvCandidates.BreakpointsFilter, sub2)
#   }
# }
# table(InvCandidates.BreakpointsFilter$fail)
# # 
# InvCandidates.BreakpointsFilter = subset(InvCandidates.BreakpointsFilter, fail == 0)[1:5]
# # write to file
# write.table(file=paste0(dir,"InvCandidates.BreakpointsFilter-2-",ref,".txt"), InvCandidates.BreakpointsFilter, quote=F, row.names=F, sep="\t")



#################### ADD INFO OF HOW CLOSE CANDIDATE INVERSIONS ####################
#################### ARE TO SCAFFOLD BORDERS WITHIN REF GENOMES ####################

## GET INFORMATION ABOUT CHROMOSOME BREAKPOINTS ==================================
# read if existent
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)

# if not read original and update
if (!exists("scaffLengths")) {
  # get scaffold start and end positions relative to chromosome ////////////////////
  scaffLengths = read.table(paste0("/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/",ref,".scaffoldsLength.txt"))
  # get chromosomes to which each ref. scaffolds belongs to 
  names(scaffLengths) = c("scaffold", "length")
  if (ref == "hmelv25") { scaffLengths$chrom = substr(scaffLengths$scaffold, 1, 7) }
  if (ref == "heradem") { scaffLengths$chrom = substr(scaffLengths$scaffold, 1, 8) }
  if (ref == "heralat") { scaffLengths$chrom = sapply(str_split(scaffLengths$scaffold, "_"), "[", 2) }
  ## update scaffold coordinates (by their order in chromosome) ====================
  scaffLengths$c_sta = 0
  scaffLengths$c_end = 0
  # define first line cumulative coordinates
  scaffLengths$c_sta[1] = 1
  scaffLengths$c_end[1] = scaffLengths$length[1]
  # update cumulative coordinates
  for (rw in 2:nrow(scaffLengths)) {
    # in the same chromosome as previous scaffold
    if (scaffLengths$chrom[rw] == scaffLengths$chrom[rw-1]) {
      scaffLengths$c_sta[rw] = scaffLengths$c_sta[rw] + scaffLengths$c_end[rw-1] + 1
      scaffLengths$c_end[rw] = scaffLengths$c_sta[rw] + scaffLengths$length[rw]
    }
    # start of a different chromosome 
    if (scaffLengths$chrom[rw] != scaffLengths$chrom[rw-1]) {
      scaffLengths$c_sta[rw] = 1
      scaffLengths$c_end[rw] = scaffLengths$length[rw]
    }
  }
  # write to file
  write.table(file=paste0(dir,ref,".scaffLengths.txt"), scaffLengths, quote=F, row.names=F)
}



## UPDATE INVERSION BREAKPOINT COORDINATES ======================================
# arrange data so that same chromosome scaffolds come together ==================
if (ref == "hmelv25") {
  InvCandidates.Breakpoints$chrom = substr(InvCandidates.Breakpoints$R.scaffold, 1, 7)
  InvCandidates.Breakpoints$scaff = substr(InvCandidates.Breakpoints$R.scaffold, 8, 11)
  InvCandidates.Breakpoints = subset(InvCandidates.Breakpoints, chrom != "Hmel200")
}
if (ref == "heradem") {
  InvCandidates.Breakpoints$chrom = substr(InvCandidates.Breakpoints$R.scaffold, 1, 8)
  InvCandidates.Breakpoints$scaff = substr(InvCandidates.Breakpoints$R.scaffold, 9, 10)
  InvCandidates.Breakpoints = subset(InvCandidates.Breakpoints, R.scaffold != "Herato_mt")
}
if (ref == "heralat") {
  InvCandidates.Breakpoints$chrom = sapply(str_split(InvCandidates.Breakpoints$R.scaffold, "_"), "[", 2)
  InvCandidates.Breakpoints$scaff = sapply(str_split(InvCandidates.Breakpoints$R.scaffold, "_"), "[", 3)
}
# update Inversion coordinates ==================================================
InvCandidates.Breakpoints$c_Bsta = 0
InvCandidates.Breakpoints$c_Bend = 0
for (rw in 1:nrow(InvCandidates.Breakpoints)) {
  scaff = as.character(InvCandidates.Breakpoints[rw,]$R.scaffold)
  csta = subset(scaffLengths, scaffold == scaff)$c_sta
  InvCandidates.Breakpoints$c_Bsta[rw] = InvCandidates.Breakpoints$Bsta[rw] + csta
  InvCandidates.Breakpoints$c_Bend[rw] = InvCandidates.Breakpoints$Bend[rw] + csta
}


## filter Inversions that are right next to scaffold boundaries =================
## [likely to be artifacts] =====================================================
InvCandidates.Breakpoints$distBorder = -1
# InvCandidates.BreakpointsFilter1 = data.frame()

for (rw in 1:nrow(InvCandidates.Breakpoints)) {
  # get corresponding scaffold left and right coordinates
  scaff = as.character(InvCandidates.Breakpoints[rw,]$R.scaffold)
  llim = subset(scaffLengths, scaffold == scaff)$c_sta
  rlim = subset(scaffLengths, scaffold == scaff)$c_end
  # get Inversion coordinates
  ista = InvCandidates.Breakpoints[rw,]$c_Bsta
  iend = InvCandidates.Breakpoints[rw,]$c_Bend
  # distance to scaffold limits
  minDist = min(abs(ista-llim), abs(ista-rlim), abs(iend-llim), abs(iend-rlim))
  InvCandidates.Breakpoints$distBorder[rw] = minDist
  # if (minDist >= minAllowedDistance) {InvCandidates.BreakpointsFilter1 = rbind(InvCandidates.BreakpointsFilter1, InvCandidates.Breakpoints[rw,])}
}
# write to file
write.table(file=paste0(dir,"InvCandidates.Breakpoints-2-",ref,".txt"), InvCandidates.Breakpoints, quote=F, row.names=F, sep="\t")



#################### ADD COVERAGE DATA #########################

# coverage data
covTable = read.table(paste0(dir,"coverageTable.",ref,".txt"), header=T)

# go through each species and apply filter
tmp = data.frame()
for (spp in unique(InvCandidates.Breakpoints$species)) {
  s1 = subset(InvCandidates.Breakpoints, species == spp)
  # go through each scaffold
  for (contig in unique(s1$Q.scaffold)) {
    s2 = subset(s1, Q.scaffold == contig)
    sf = as.character(s2$R.scaffold)
    st = s2$Bsta
    en = s2$Bend
    stCov = subset(covTable, species == spp & scaffold == sf & sta < st & end > st)$relcov
    enCov = subset(covTable, species == spp & scaffold == sf & sta < en & end > en)$relcov
    if (length(stCov) > 0) {s2$stCov = stCov}; if (length(stCov) == 0) {s2$stCov = NA} 
    if (length(enCov) > 0) {s2$enCov = enCov}; if (length(enCov) == 0) {s2$enCov = NA}
    tmp = rbind(tmp, s2)
  }
}
InvCandidates.Breakpoints = tmp
# write to file
write.table(file=paste0(dir,"InvCandidates.Breakpoints-2-",ref,".txt"), InvCandidates.Breakpoints, quote=F, row.names=F)

# plot coverage
ggplot(tmp) + 
  geom_boxplot(aes(x=species, y=stCov)) +
  geom_boxplot(aes(x=species, y=enCov)) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1))



######################## APPLY LAST FILTER ##############################
## read original minimap2 mappings =====
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
subdir = "2.minimap2"
filenames <- list.files(paste0(dir, subdir), pattern=paste0("*",ref,".filter.mm2"), full.names=TRUE)
# read data
if (!exists("minimaps")) {
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
}

## filter out candidate breakpoints that have a "straight" alignment going through it =====
InvCandidates.Breakpoints2 = data.frame()
for (rw in 1:nrow(InvCandidates.Breakpoints)) {
  #
  state = "pass"; ov1 = 0; ov2 = 0
  #
  s1 = InvCandidates.Breakpoints[rw,]
  # define variables
  spp = s1$species
  rsc = s1$R.scaffold
  bk1 = s1$Bsta
  bk2 = s1$Bend
  # get mappings of [spp] to the [refScaffold]
  s2 = subset(minimaps, species == spp & R.scaffold == as.character(rsc) & MQ >= 10)
  s2.bk1 = subset(s2, R.start < bk1 & R.end > bk1); 
  s2.bk2 = subset(s2, R.start < bk2 & R.end > bk2);
  if (nrow(s2.bk1) > 0) {
    s2.bk1$overlap = -1
    for (i in nrow(s2.bk1)) {
      s2.bk1$overlap[i] = min(abs(bk1 - s2.bk1$R.start[i]), abs(bk1 - s2.bk1$R.end[i]))
    }
  }
  if (nrow(s2.bk2) > 0) {
    s2.bk2$overlap = -1
    for (i in nrow(s2.bk2)) {
      s2.bk2$overlap[i] = min(abs(bk2 - s2.bk2$R.start[i]), abs(bk2 - s2.bk2$R.end[i]))
    }
  }
  s2.cmb = rbind(s2.bk1, s2.bk2)
  if (nrow(s2.cmb) > 0) {
    s2.cmbOverlap = subset(s2.cmb, overlap >= 500)
    if (nrow(s2.cmbOverlap) >= 1) {state = "fail"}
    #if (ov1 >= 500 | ov2 >= 500) {state = "fail"}
    # if (ov1+ov2 >= 1000) {state = "fail"}
  }
  if (state == "pass") { InvCandidates.Breakpoints2 = rbind(InvCandidates.Breakpoints2, s1) }
}
# write to file
write.table(file=paste0(dir,"InvCandidates.Breakpoints2-2-",ref,".txt"), InvCandidates.Breakpoints2, quote=F, row.names=F)



# chr = bigchroms[17]
# spSummaryFilterSub = subset(InvCandidates.BreakpointsFilter4, chrom == chr)
# 
# scaffLengthsSub = subset(scaffLengths, chrom == chr)
# # plot all species
# ggplot(spSummaryFilterSub) +
#   geom_tile(aes(x=(c_Bsta+c_Bend)/2000000, y=species, width=(c_Bend-c_Bsta)/1000000, height=0.8, fill=species), col="black", alpha=0.5) +
#   geom_vline(data=scaffLengthsSub, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
#   geom_vline(data=scaffLengthsSub, aes(xintercept = c_end/1000000), col="grey", lty=2) +
#   geom_text(aes(x=c_Bsta/1000000, y=species, label=round(stCov,2))) +
#   geom_text(aes(x=c_Bend/1000000, y=species, label=round(enCov,2))) +
#   # xlim(xmin,xmax) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
#   ggtitle(chr) +
#   xlab("Chromosome Position(Mb)") + ylab("Species")


# #################### Filter candidates based on size ####################
# minInvSize = 50000
# 
# # filter candidates based on size
# # InvCandidates.BreakpointsFilter2 = subset(InvCandidates.BreakpointsFilter, Bend-Bsta >= minInvSize)
# InvCandidates.BreakpointsFilter2 = subset(InvCandidates.Breakpoints, Bend-Bsta >= minInvSize)
# InvCandidates.BreakpointsFilter2 = subset(InvCandidates.BreakpointsFilter2, Bend-Bsta < 3000000)
# # plot size of Inversion Candidates
# ggplot(InvCandidates.BreakpointsFilter2, aes(x=species, y=(Bend-Bsta)/1000)) +
#   geom_boxplot(alpha=0.25) + geom_jitter(alpha=0.5) + coord_flip() +
#   xlab("Species") + ylab("Candidate Inversion Size (Kb)")
# # write to file
# write.table(file=paste0(dir,"InvCandidates.BreakpointsFilter2-2-",ref,".txt"), InvCandidates.BreakpointsFilter2, quote=F, row.names=F, sep="\t")
# 
# #################### Output InvCandidates Query coordinates ####################
# # create dataframe
# InvCand.QueryCoords = data.frame()
# # parse
# for (spp in unique(InvCandidates.Merge2$species)) {
#   sub1 = subset(InvCandidates.Merge2, species == spp)
#   for (contig in unique(sub1$Q.scaffold)) {
#     sub2 = subset(sub1, Q.scaffold == contig)
#     pl = subset(sub2, strand == "+")
#     mn = subset(sub2, strand == "-")
#     df = data.frame(species=spp, Scaffold=contig, sta1=pl$Q.start, end1=pl$Q.end, sta2=mn$Q.start, end2=mn$Q.end)
#     InvCand.QueryCoords = rbind(InvCand.QueryCoords, df)
#   }
# }
# # write to file
# write.table(file=paste0(dir,"InvCand.QueryCoords-2-",ref,".txt"), InvCand.QueryCoords, quote=F, row.names=F, sep="\t")
