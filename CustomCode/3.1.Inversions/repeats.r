
## read data ==================================================
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
fil = "hmelv25.repeatInfo.txt"
Repeats = read.table(paste0(dir,fil))
names(Repeats) = c("scaffold", "start", "end", "RepeatType", "percDiv", "PercDel", "PercIns")

# add chromosome names
if (ref == "hmelv25") {
  Repeats$chrom = substr(Repeats$scaffold, 1, 7)
  Repeats = subset(Repeats, chrom != "Hmel200")
}
if (ref == "heradem") {
  Repeats$chrom = substr(Repeats$R.scaffold, 1, 8)
  Repeats = subset(Repeats, R.scaffold != "Herato_mt")
}



## update repeat coordinates ==================================
# ref chromosome coordinates
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)

# subset to chromosomes
if (ref == "hmelv25") { Repeats = subset(Repeats, chrom %in% paste0("Hmel2", sprintf("%02d", seq(1,21)))) }
if (ref == "heradem") { Repeats = subset(Repeats, chrom %in% paste0("Herato", sprintf("%02d", seq(1,21)))) }
if (ref == "heralat") { Repeats = subset(Repeats, chrom %in% paste0("chr", seq(1,21))) }
# update
Repeats$c_Bsta = 0
Repeats$c_Bend = 0
for (scaff in unique(scaffLengths$scaffold)) {
  s1 = subset(scaffLengths, scaffold == scaff)
  Repeats$c_Bsta = ifelse(Repeats$scaffold == as.character(scaff), Repeats$sta + s1$c_sta, Repeats$c_Bsta)
  Repeats$c_Bend = ifelse(Repeats$scaffold == as.character(scaff), Repeats$end + s1$c_sta - 1, Repeats$c_Bend)
}
