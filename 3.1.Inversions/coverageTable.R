## load modules
require(stringr)

## variables
ref = "heradem"
dir=paste0("/n/mallet_lab/fseixas/1.projects/1.pseudo_references/3.1.compareMap/1.2.rmdup/",ref,"/")
files = list.files(path = dir, patter = ".cov.")

## get coverage data
covTable = data.frame()
for (file in files) {
  inp = read.table(paste0(dir,file))
  names(inp) = c("scaffold","sta","end","sites","cov","reference")
  inp$relcov = inp$cov/median(inp$cov)
  inp$pos = seq(1,nrow(inp),1)
  inp$relpos = inp$pos/max(inp$pos)
  covTable = rbind(covTable, inp)
}

## update species names
covTable$reference = tolower(covTable$reference)
covTable$reference = str_remove(covTable$reference, "_")
covTable$reference = str_remove(covTable$reference, "_")
covTable$reference = str_replace(covTable$reference, "hel", "h")
covTable$reference = str_replace(covTable$reference, "lapdor", "hdor")
unique(covTable$reference)
names(covTable)[6] = "species"

## define chromosomes
if (ref == "hmelv25") {covTable$chrom = substr(covTable$scaffold, 1, 7)}
if (ref == "heradem") {covTable$chrom = substr(covTable$scaffold, 1, 8)}
# if (ref == "heralat") {covTable$chrom = substr(covTable$scaffold, 1, 7)}


# update coverage Table coordinates ==================================================
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
# ref chromosome coordinates
scaffLengths = read.table(paste0(dir,ref,".scaffLengths.txt"), header=T)
# subset to chromosomes
if (ref == "hmelv25") { covTable = subset(covTable, chrom %in% paste0("Hmel2", sprintf("%02d", seq(1,21)))) }
if (ref == "heradem") { covTable = subset(covTable, chrom %in% paste0("Herato", sprintf("%02d", seq(1,21)))) }
if (ref == "heralat") { covTable = subset(covTable, chrom %in% paste0("chr", seq(1,21))) }
# update
covTable$c_Bsta = 0
covTable$c_Bend = 0
for (scaff in unique(scaffLengths$scaffold)) {
  s1 = subset(scaffLengths, scaffold == scaff)
  covTable$c_Bsta = ifelse(covTable$scaffold == as.character(scaff), covTable$sta + s1$c_sta, covTable$c_Bsta)
  covTable$c_Bend = ifelse(covTable$scaffold == as.character(scaff), covTable$end + s1$c_sta - 1, covTable$c_Bend)
}

## write to file
dir="/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/"
write.table(file=paste0(dir, "coverageTable.",ref,".txt"), covTable, row.names = F, quote=F)


