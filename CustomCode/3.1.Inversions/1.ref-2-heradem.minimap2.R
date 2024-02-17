## load modules
require(stringr)


## variables
ref="heradem"
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
# ref chromosome coordinates
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)


#################### hmelv25-2-heradem ####################
hmelv25_2_heradem_Minimap2 = read.table(paste0(dir,"2.minimap2/hmelv25-2-",ref,".filter.mm2"))
names(hmelv25_2_heradem_Minimap2)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
# stick to chromosomes
if (ref == "heradem") {
  hmelv25_2_heradem_Minimap2$chrom = substr(hmelv25_2_heradem_Minimap2$R.scaffold, 1, 8)
  hmelv25_2_heradem_Minimap2 = subset(hmelv25_2_heradem_Minimap2, Q.scaffold != "Herato_mt")
}
# update hmelv25-2-heradem coordinates ===========
hmelv25_2_heradem_Minimap2$c_Bsta = 0
hmelv25_2_heradem_Minimap2$c_Bend = 0
for (scaff in unique(scaffLengths$scaffold)) {
  csta = subset(scaffLengths, scaffold == scaff)$c_sta
  hmelv25_2_heradem_Minimap2$c_Bsta = ifelse(hmelv25_2_heradem_Minimap2$R.scaffold == scaff, hmelv25_2_heradem_Minimap2$R.start + csta, hmelv25_2_heradem_Minimap2$c_Bsta)
  hmelv25_2_heradem_Minimap2$c_Bend = ifelse(hmelv25_2_heradem_Minimap2$R.scaffold == scaff, hmelv25_2_heradem_Minimap2$R.end + csta, hmelv25_2_heradem_Minimap2$c_Bend)
}
# filter by MQ and matching chromosome
hmelv25_2_heradem_Minimap2$Qchr = substr(hmelv25_2_heradem_Minimap2$Q.scaffold, 6, 7)
hmelv25_2_heradem_Minimap2$Rchr = substr(hmelv25_2_heradem_Minimap2$R.scaffold, 7, 8)
hmelv25_2_heradem_Minimap2 = subset(hmelv25_2_heradem_Minimap2, MQ >= 20 & Qchr == Rchr)

# merge alignments
hmelv25_2_heradem_Minimap2Merge = mergeAligns(hmelv25_2_heradem_Minimap2)

# update R coordinates
hmelv25_2_heradem_Minimap2Merge = updateRcoords(hmelv25_2_heradem_Minimap2Merge)

# write to file
write.table(file=paste0(dir,"hmelv25_2_heradem_Minimap2.txt"), hmelv25_2_heradem_Minimap2, quote=F, row.names=F)
write.table(file=paste0(dir,"hmelv25_2_heradem_Minimap2Merge.txt"), hmelv25_2_heradem_Minimap2Merge, quote=F, row.names=F)




#################### HERALAT-2-heradem ####################
heralat_2_heradem_Minimap2 = read.table(paste0(dir,"2.minimap2/heralat-2-",ref,".filter.mm2"))
names(heralat_2_heradem_Minimap2)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
head(heralat_2_heradem_Minimap2)

# stick to chromosomes
if (ref == "heradem") {
  heralat_2_heradem_Minimap2$chrom = substr(heralat_2_heradem_Minimap2$R.scaffold, 1, 8)
  heralat_2_heradem_Minimap2 = subset(heralat_2_heradem_Minimap2, Q.scaffold != "Herato_mt")
}

# update heralat-2-heradem coordinates ===========
heralat_2_heradem_Minimap2$c_Bsta = 0
heralat_2_heradem_Minimap2$c_Bend = 0
for (scaff in unique(scaffLengths$scaffold)) {
  csta = subset(scaffLengths, scaffold == scaff)$c_sta
  heralat_2_heradem_Minimap2$c_Bsta = ifelse(heralat_2_heradem_Minimap2$R.scaffold == scaff, heralat_2_heradem_Minimap2$R.start + csta, heralat_2_heradem_Minimap2$c_Bsta)
  heralat_2_heradem_Minimap2$c_Bend = ifelse(heralat_2_heradem_Minimap2$R.scaffold == scaff, heralat_2_heradem_Minimap2$R.end + csta, heralat_2_heradem_Minimap2$c_Bend)
}

# filter by MQ and matching chromosome
heralat_2_heradem_Minimap2$Qchr = sapply(str_split(heralat_2_heradem_Minimap2$Q.scaffold, "_"), "[", 2)
heralat_2_heradem_Minimap2$Qchr = str_remove(heralat_2_heradem_Minimap2$Qchr, "chr")
heralat_2_heradem_Minimap2$Qchr = ifelse(heralat_2_heradem_Minimap2$Qchr %in% seq(1,9), paste0("0", heralat_2_heradem_Minimap2$Qchr), heralat_2_heradem_Minimap2$Qchr)
heralat_2_heradem_Minimap2$Rchr = substr(heralat_2_heradem_Minimap2$R.scaffold, 7, 8)
heralat_2_heradem_Minimap2 = subset(heralat_2_heradem_Minimap2, MQ >= 20 & Qchr == Rchr)

# merge alignments
heralat_2_heradem_Minimap2Merge = mergeAligns(heralat_2_heradem_Minimap2)

# update R.coordinates and C.coordinates
heralat_2_heradem_Minimap2Merge = updateRcoords(heralat_2_heradem_Minimap2Merge)

# write to file
write.table(file=paste0(dir,"heralat_2_heradem_Minimap2.txt"), heralat_2_heradem_Minimap2, quote=F, row.names=F)
write.table(file=paste0(dir,"heralat_2_heradem_Minimap2Merge.txt"), heralat_2_heradem_Minimap2Merge, quote=F, row.names=F)

nrow(hmelv25_2_heradem_Minimap2Merge)
nrow(heralat_2_heradem_Minimap2Merge)




# aaa = subset(hmelv25_2_heradem_Minimap2Merge, chrom == "Herato06")
# aaa = subset(heralat_2_heradem_Minimap2Merge, chrom == "Herato06")
# tail(aaa)
# ggplot(aaa) +
#   # geom_gene_arrow(data=aaa, aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=0), fill=ifelse(aaa$strand == "-", "#b2182b", "#f5f5f5"), lwd=0.6) +
#   geom_gene_arrow(data=aaa, aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=0), fill=ifelse(aaa$strand == "-", "#b2182b", "#f5f5f5"), lwd=0.6) +
#   # add ref scaffold limits
#   geom_vline(data=sub_scaffLen, aes(xintercept = c_sta/1000000), col="grey", lty=2) +
#   geom_vline(data=sub_scaffLen, aes(xintercept = c_end/1000000), col="grey", lty=2) +
#   facet_wrap(~R.scaffold, nrow=7)
# 
