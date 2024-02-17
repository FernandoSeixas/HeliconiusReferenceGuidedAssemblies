## read files ==================================================
ref = "heradem"
if (ref == "heradem") {dir = "/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hera.dem.v1/"}
if (ref == "hmelv25") {dir = "/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmel2.5/"}
group = paste0(ref, ".scaffGroups")
sizes = paste0(ref,  ".scaffold.sizes")
genes = paste0(dir, ref, ".genes_coord.tab")

grp = read.table(group); names(grp) = c("Chrom", "Scaffold")
siz = read.table(sizes); names(siz) = c("Scaffold", "Size")
gen = read.table(genes); names(gen) = c("gene","Scaffold","sta","end","strand")


## determine new start coordinates of scaffold =================
chr2scaSizes = merge(grp, siz, id.vars="Scaffold")
chr2scaSizes$newStart = 1
preChrom = chr2scaSizes$Chrom[1]
preSize = chr2scaSizes$Size[1]
preCumSize = chr2scaSizes$Size[1]

for (rw in 2:nrow(chr2scaSizes)) {
  # in same chromosome
  if (chr2scaSizes$Chrom[rw] == preChrom) {
    chr2scaSizes$newStart[rw] = preCumSize + 101
    preCumSize = preCumSize + chr2scaSizes$Size[rw] + 100
  }
  # in diff chromosome
  if (chr2scaSizes$Chrom[rw] != preChrom) {
    preChrom = chr2scaSizes$Chrom[rw]
    preCumSize = chr2scaSizes$Size[rw]
  }
}
# chr2scaSizes[16,3]+chr2scaSizes[16,4]


## update gene coordinates =====================================
gen$Scaffold = as.character(gen$Scaffold)
for (rw in 1:nrow(gen)) {
  sca = as.character(gen$Scaffold[rw])
  sub = subset(chr2scaSizes, Scaffold == sca)
  nchr = as.character(sub$Chrom)
  plus = sub$newStart
  gen$Scaffold[rw] = nchr
  gen$sta[rw] = gen$sta[rw] + plus - 1
  gen$end[rw] = gen$end[rw] + plus - 1
}
write.table(file=paste0(ref,".genes_coord.vGlued.tab"), gen, quote=F, row.names = F, col.names = F, sep="\t")
