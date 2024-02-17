## load libraries
require(stringr)
require(ggrepel)

## Glued genomes - big chromosomes list
hmelv25BigChroms = as.character(read.table("/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/hmelv25.vGlued.scaffolds")[1:21,1])
herademBigChroms = as.character(read.table("/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/heradem.vGlued.scaffolds")[1:21,1])
heralatBigChroms = as.character(read.table("/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/0.1.references/heralat.vGlued.scaffolds")[1:21,1])

## 
# read scaffold lengths table
hmelv25ScaLens = read.table("/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmel2.5/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.lengths.txt")
herademScaLens = read.table("/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hera.dem.v1/Heliconius_erato_demophoon_v1_-_scaffolds.lengths.txt") 
heralatScaLens = read.table("/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hera.lat.v1/Heliconius_erato_lativitta_v1_-_scaffolds.lengths.txt") 
names(hmelv25ScaLens) = c("Scaffold","Length")
names(herademScaLens) = c("Scaffold","Length")
names(heralatScaLens) = c("Scaffold","Length")
# separate chromosome and scaffold number
hmelv25ScaLens$chrom = substr(hmelv25ScaLens$Scaffold, 1, 7)
herademScaLens$chrom = substr(herademScaLens$Scaffold, 1, 8)
heralatScaLens$chrom = str_replace(str_remove_all(substr(heralatScaLens$Scaffold, 1, 9),"_"),"Hel","Hel_")
hmelv25ScaLens$ScaNb = substr(hmelv25ScaLens$Scaffold, 8, 10)
herademScaLens$ScaNb = substr(herademScaLens$Scaffold, 9, 11)
heralatScaLens$ScaNb = str_remove(substr(heralatScaLens$Scaffold, 10, 11),"_")
# restrict to big chromosomes
hmelv25ScaLens = subset(hmelv25ScaLens, chrom %in% hmelv25BigChroms)
herademScaLens = subset(herademScaLens, chrom %in% herademBigChroms)
heralatScaLens = subset(heralatScaLens, chrom %in% heralatBigChroms)

## cumulative lengths ==================================================
# hmelv25
df = data.frame()
for (ch in unique(hmelv25ScaLens$chrom)) {
  s1 = subset(hmelv25ScaLens, chrom == ch)
  s1$cumSta = 1
  s1$cumEnd = 0
  s1$cumEnd[1] = s1$Length[1]
  if (nrow(s1) > 1) {
    for (rw in 2:nrow(s1)) {
      s1$cumSta[rw] = s1$cumEnd[(rw-1)]+101
      s1$cumEnd[rw] = s1$cumSta[rw] + s1$Length[rw] - 1
    }
  }
  df = rbind(df, s1)
}
hmelv25ScaLens = df

# heradem
df = data.frame()
for (ch in unique(herademScaLens$chrom)) {
  s1 = subset(herademScaLens, chrom == ch)
  s1$cumSta = 1
  s1$cumEnd = 0
  s1$cumEnd[1] = s1$Length[1]
  if (nrow(s1) > 1) {
    for (rw in 2:nrow(s1)) {
      s1$cumSta[rw] = s1$cumEnd[(rw-1)]+101
      s1$cumEnd[rw] = s1$cumSta[rw] + s1$Length[rw] - 1
    }
  }
  df = rbind(df, s1)
}
herademScaLens = df

# heralat
df = data.frame()
for (ch in unique(heralatScaLens$chrom)) {
  s1 = subset(heralatScaLens, chrom == ch)
  s1$cumSta = 1
  s1$cumEnd = 0
  s1$cumEnd[1] = s1$Length[1]
  if (nrow(s1) > 1) {
    for (rw in 2:nrow(s1)) {
      s1$cumSta[rw] = s1$cumEnd[(rw-1)]+101
      s1$cumEnd[rw] = s1$cumSta[rw] + s1$Length[rw] - 1
    }
  }
  df = rbind(df, s1)
}
heralatScaLens = df

write.table(file=paste0(dir, "hmelv25.ScaffoldBreakpoints.txt"), hmelv25ScaLens, quote=F, row.names=F)
write.table(file=paste0(dir, "heradem.ScaffoldBreakpoints.txt"), herademScaLens, quote=F, row.names=F)
write.table(file=paste0(dir, "heralat.ScaffoldBreakpoints.txt"), heralatScaLens, quote=F, row.names=F)


## convert coordinates
rf = "hmelv25"
ch = "Hmel207"
st = 5852090
en = 5916429
en-st
if (rf == "hmelv25") { s1 = subset(hmelv25ScaLens, chrom == ch) }
if (rf == "heradem") { s1 = subset(herademScaLens, chrom == ch) }
if (rf == "heralat") { s1 = subset(heralatScaLens, chrom == ch) }

s1 = subset(s1, en > cumSta & st < cumEnd)

if (nrow(s1) == 1) {
  c = s1$Scaffold
  s = st-s1$cumSta+1
  e = en-s1$cumSta+1
  print(paste0(c,":",s,"-",e))
}

if (nrow(s1) > 1) {
  s = subset(s1, st > cumSta & st < cumEnd)
  e = subset(s1, en > cumSta & en < cumEnd)
  cl = s$Scaffold
  sl = s$Length-(s$cumEnd-st);
  el = s$Length
  cr = e$Scaffold
  sr = 1
  er = en-e$cumSta
  print(s1)
  print(paste0(cl,":",sl,"-",el))
  print(paste0(cr,":",sr,"-",er))
}


## PCA
dir1 = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.1.Inversions/minimap2/Inv2IQtree/"
vec = read.table(paste0(dir1,"Hmel219002o.1.147605.hmelv25_clade_all.eigenvec"))
names(vec) = c("POP","IND",paste0("PC",seq(1,ncol(vec)-3)))

# clade assignment
mel_clade = c("hmel","hcyd","htim","hbes","hnum","hhec","hele","hpar")
mid_clade = c("hbur","hdor","haoe")
era_clade = c("hera","hhimfat","hhim","hsia","htel","hdem","hsar")

vec$IND = tolower(str_replace(vec$IND, "HEL_", "h"))
vec$IND = substr(vec$IND, 1, 4)
vec$IND = ifelse(vec$IND == "lap_", "hdor", vec$IND)
vec$IND = ifelse(vec$IND == "ner_", "haoe", vec$IND)

vec$POP = ifelse(vec$IND %in% mel_clade, "mel", as.character(vec$POP))
vec$POP = ifelse(vec$IND %in% mid_clade, "mid", as.character(vec$POP))
vec$POP = ifelse(vec$IND %in% era_clade, "era", as.character(vec$POP))

# plot
ggplot(vec, aes(x=PC1, y=PC2, colour=POP, label=IND)) + geom_point() + geom_text_repel()

