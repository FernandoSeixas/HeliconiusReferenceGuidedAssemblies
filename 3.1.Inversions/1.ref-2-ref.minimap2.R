## load modules
require(stringr)


## variables
ref="hmelv25"
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
# ref chromosome coordinates
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)


#################### HERADEM-2-HMELV25 ####################
heradem_2_hmelv25_Minimap2 = read.table(paste0(dir,"2.minimap2/heradem-2-",ref,".filter.mm2"))
names(heradem_2_hmelv25_Minimap2)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
# stick to chromosomes
if (ref == "hmelv25") {
  heradem_2_hmelv25_Minimap2$chrom = substr(heradem_2_hmelv25_Minimap2$R.scaffold, 1, 7)
  heradem_2_hmelv25_Minimap2 = subset(heradem_2_hmelv25_Minimap2, chrom != "Hmel200")
  heradem_2_hmelv25_Minimap2 = subset(heradem_2_hmelv25_Minimap2, Q.scaffold != "Herato_mt")
}
# update heradem-2-hmelv25 coordinates ===========
heradem_2_hmelv25_Minimap2$c_Bsta = 0
heradem_2_hmelv25_Minimap2$c_Bend = 0
for (scaff in unique(scaffLengths$scaffold)) {
  csta = subset(scaffLengths, scaffold == scaff)$c_sta
  heradem_2_hmelv25_Minimap2$c_Bsta = ifelse(heradem_2_hmelv25_Minimap2$R.scaffold == scaff, heradem_2_hmelv25_Minimap2$R.start + csta, heradem_2_hmelv25_Minimap2$c_Bsta)
  heradem_2_hmelv25_Minimap2$c_Bend = ifelse(heradem_2_hmelv25_Minimap2$R.scaffold == scaff, heradem_2_hmelv25_Minimap2$R.end + csta, heradem_2_hmelv25_Minimap2$c_Bend)
}
# filter by MQ and matching chromosome
heradem_2_hmelv25_Minimap2$Qchr = substr(heradem_2_hmelv25_Minimap2$Q.scaffold, 7, 8)
heradem_2_hmelv25_Minimap2$Rchr = substr(heradem_2_hmelv25_Minimap2$R.scaffold, 6, 7)
heradem_2_hmelv25_Minimap2 = subset(heradem_2_hmelv25_Minimap2, MQ >= 20 & Qchr == Rchr)

# merge alignments
heradem_2_hmelv25_Minimap2Merge = mergeAligns(heradem_2_hmelv25_Minimap2)

# update R coordinates
heradem_2_hmelv25_Minimap2Merge = updateRcoords(heradem_2_hmelv25_Minimap2Merge)

# write to file
write.table(file=paste0(dir,"heradem_2_hmelv25_Minimap2.txt"), heradem_2_hmelv25_Minimap2, quote=F, row.names=F)
write.table(file=paste0(dir,"heradem_2_hmelv25_Minimap2Merge.txt"), heradem_2_hmelv25_Minimap2Merge, quote=F, row.names=F)



#################### HERALAT-2-HMELV25 ####################
heralat_2_hmelv25_Minimap2 = read.table(paste0(dir,"2.minimap2/heralat-2-",ref,".filter.mm2"))
names(heralat_2_hmelv25_Minimap2)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
head(heralat_2_hmelv25_Minimap2)
# stick to chromosomes
if (ref == "hmelv25") {
  heralat_2_hmelv25_Minimap2$chrom = substr(heralat_2_hmelv25_Minimap2$R.scaffold, 1, 7)
  heralat_2_hmelv25_Minimap2 = subset(heralat_2_hmelv25_Minimap2, chrom != "Hmel200")
}
head(heralat_2_hmelv25_Minimap2)

# update heralat-2-hmelv25 coordinates ===========
heralat_2_hmelv25_Minimap2$c_Bsta = 0
heralat_2_hmelv25_Minimap2$c_Bend = 0
for (scaff in unique(scaffLengths$scaffold)) {
  csta = subset(scaffLengths, scaffold == scaff)$c_sta
  heralat_2_hmelv25_Minimap2$c_Bsta = ifelse(heralat_2_hmelv25_Minimap2$R.scaffold == scaff, heralat_2_hmelv25_Minimap2$R.start + csta, heralat_2_hmelv25_Minimap2$c_Bsta)
  heralat_2_hmelv25_Minimap2$c_Bend = ifelse(heralat_2_hmelv25_Minimap2$R.scaffold == scaff, heralat_2_hmelv25_Minimap2$R.end + csta, heralat_2_hmelv25_Minimap2$c_Bend)
}
# filter by MQ and matching chromosome
heralat_2_hmelv25_Minimap2$Qchr = sapply(str_split(heralat_2_hmelv25_Minimap2$Q.scaffold, "_"), "[", 2)
heralat_2_hmelv25_Minimap2$Qchr = str_remove(heralat_2_hmelv25_Minimap2$Qchr, "chr")
heralat_2_hmelv25_Minimap2$Qchr = ifelse(heralat_2_hmelv25_Minimap2$Qchr %in% seq(1,9), paste0("0", heralat_2_hmelv25_Minimap2$Qchr), heralat_2_hmelv25_Minimap2$Qchr)
heralat_2_hmelv25_Minimap2$Rchr = substr(heralat_2_hmelv25_Minimap2$R.scaffold, 6, 7)
heralat_2_hmelv25_Minimap2 = subset(heralat_2_hmelv25_Minimap2, MQ >= 20 & Qchr == Rchr)

# merge alignments
heralat_2_hmelv25_Minimap2Merge = mergeAligns(heralat_2_hmelv25_Minimap2)

# update R coordinates
heralat_2_hmelv25_Minimap2Merge = updateRcoords(heralat_2_hmelv25_Minimap2Merge)

# write to file
write.table(file=paste0(dir,"heralat_2_hmelv25_Minimap2.txt"), heralat_2_hmelv25_Minimap2, quote=F, row.names=F)
write.table(file=paste0(dir,"heralat_2_hmelv25_Minimap2Merge.txt"), heralat_2_hmelv25_Minimap2Merge, quote=F, row.names=F)

nrow(heradem_2_hmelv25_Minimap2Merge)
nrow(heralat_2_hmelv25_Minimap2Merge)
