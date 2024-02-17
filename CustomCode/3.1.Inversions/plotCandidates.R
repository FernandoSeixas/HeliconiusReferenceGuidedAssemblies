# load libraries
require(ggplot2)

# variables
ref="heradem"
dir="/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/vGlued/"
# define clades
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")
# reference assembly scaffold breakpoints
refBreaks = read.table(paste0(dir, "../", ref, ".ScaffoldBreakpoints.txt"), header=T)

#
InvCandidates.BreakpointsFilter = read.table(paste0(dir,"InvCandidates.BreakpointsFilter-2-",ref,".txt"), header=T)
InvCandidates.BreakpointsFilter$species = factor(InvCandidates.BreakpointsFilter$species, levels = c(mel_clade,mid_clade,era_clade))

ppp1 = subset(InvCandidates.BreakpointsFilter, Bend-Bsta >= 50000)
ppp2 = subset(ppp1, Bend-Bsta < 3000000)
# table(InvCandidates.BreakpointsFilter$species)
# ggplot(ppp2, aes(x=species, y=(Bend-Bsta)/1000)) +
#   geom_boxplot(alpha=0.25) +
#   geom_jitter(alpha=0.5) +
#   coord_flip()

##
df1 = data.frame()
for (spp in unique(ppp2$species)) {
  s1 = subset(ppp2, species == spp)
  for (ch in unique(s1$R.scaffold)) {
    s2 = subset(s1, R.scaffold == ch)
    s2 = s2[order(s2$Bsta),]
    s2$count = 1
    s3 = s2[1,]
    if (nrow(s2) >= 2) { 
      for (rw in 2:nrow(s2)) {
        if (s2$Bsta[rw]-s3$Bsta[nrow(s3)] < 2500 & s2$Bend[rw]-s3$Bend[nrow(s3)] < 2500) {
          s3$Bsta[nrow(s3)] = min(s2$Bsta[rw],s3$Bsta[nrow(s3)])
          s3$Bend[nrow(s3)] = max(s2$Bend[rw],s3$Bend[nrow(s3)])
          s3$count = s3$count + 1
        }
        if (s2$Bsta[rw]-s3$Bsta[nrow(s3)] >= 2500 | s2$Bend[rw]-s3$Bend[nrow(s3)] >= 2500) {
          s3 = rbind(s3, s2[rw,])
        }
      } 
    }
    df1 = rbind(df1, s3)
  }
}
df1$len = df1$Bend - df1$Bsta + 1

# table(df1$R.scaffold)
# table(subset(df1, count >= 2)$R.scaffold)
# table(subset(df1, count >= 2)$R.scaffold)
# table(df1$R.scaffold)

## remove candidates overlapping more than 1 reference scaffold ========================================
df2 = data.frame()
for (rw in 1:nrow(df1)) {
  s1 = df1[rw,]
  refsub = subset(refBreaks, chrom == as.character(s1$R.scaffold))
  state = 0
  for (scn in 1:nrow(refsub)) {
    refsub2 = refsub[scn,]
    if (refsub2$cumSta > s1$Bsta & refsub2$cumSta < s1$Bend) {state = 1}
    if (refsub2$cumEnd > s1$Bsta & refsub2$cumEnd < s1$Bend) {state = 1}
  }
  if (state == 0) {df2 = rbind(df2, s1)}
}


# plot all chromosomes 
ggplot(df2, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  geom_vline(data=refBreaksSub, xintercept = refBreaksSub$cumSta/1000000) +
  # scale_x_continuous(limits=c(0,xmax), breaks=seq(0,xmax,0.50)) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggtitle(chr) + xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~R.scaffold, nrow=6)
# plot specific chromosome
chr = "Herato20"
df3 = subset(df2, R.scaffold == chr)
refBreaksSub = subset(refBreaks, chrom == chr)
ggplot(df3, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  geom_vline(data=refBreaksSub, xintercept = refBreaksSub$cumSta/1000000) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggtitle(chr) + xlab("Chromosome Position(Mb)") + ylab("Species")
df3[order(df3$Bsta),]




## group inversions concatenating all species evidence ==================================================
spSummary = data.frame()
buffer = 50000
for (ch in unique(df2$R.scaffold)) {
  s1 = subset(df2, R.scaffold == ch)
  s1 = s1[order(s1$Bsta),]
  s1$count = 0
  if (nrow(s1) >= 1) { 
    for (rw in 1:nrow(s1)) {
      st = s1$Bsta[rw]
      en = s1$Bend[rw]
      ss = subset(s1, Bsta > st-buffer & Bsta < st+buffer)
      ss = subset(ss, Bend > en-buffer & Bend < en+buffer)
      s1$count[rw] = nrow(ss)
    }
  }
  spSummary = rbind(spSummary, s1)
}
# spSummaryFilter = subset(spSummary, count > 1 & count <= 8)
spSummaryFilter = subset(spSummary, count >= 3)
table(spSummaryFilter$R.scaffold)

# plot all chromosomes
ggplot(spSummaryFilter, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  xlab("Chromosome Position(Mb)") + ylab("Species") +
  facet_wrap(~R.scaffold, nrow=5)

# plot specific chromosome
chr = "Herato11"
spSummaryFilterSub = subset(spSummaryFilter, R.scaffold == chr)

# plot all species
xmin = (min(spSummaryFilterSub$Bsta)-50000)/1000000
xmax = (max(spSummaryFilterSub$Bend)+50000)/1000000
ggplot(spSummaryFilterSub, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  geom_vline(data=refBreaksSub, xintercept = refBreaksSub$cumSta/1000000) +
  xlim(xmin,xmax) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggtitle(chr) + xlab("Chromosome Position(Mb)") + ylab("Species")

spSummaryFilterSub[order(spSummaryFilterSub$Bsta),]

# plot subset of species
sppFilter = subset(spSummaryFilterSub, species %in% c("hsar","hdem","htel","hsia","hhim","hhimfat","hera"))
# sppFilter = subset(spSummaryFilterSub, species %in% c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar"))
xmin = (min(sppFilter$Bsta)-50000)/1000000
xmax = (max(sppFilter$Bend)+50000)/1000000

ggplot(sppFilter, aes(x = (Bsta+Bend)/2000000, y=species, width=(Bend-Bsta)/1000000, height=0.8, fill=species, label=Q.scaffold)) +
  geom_tile(col="black", alpha=0.5) +
  geom_vline(data=refBreaksSub, xintercept = refBreaksSub$cumSta/1000000) +
  xlim(xmin, xmax) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggtitle(chr) + xlab("Chromosome Position(Mb)") + ylab("Species")

sppFilter[order(sppFilter$Bsta),]
