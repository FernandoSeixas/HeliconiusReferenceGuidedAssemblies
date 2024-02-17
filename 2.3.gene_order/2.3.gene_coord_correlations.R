# load libraries ----
require(optparse)

## parameters ----
option_list <- list(
  make_option(c('-r','--reference'),  action='store', default=NULL, help='reference'),
  make_option(c('-s','--species'),    action='store', default=NULL, help='species'),
  make_option(c('-g','--genecoords'),  action='store', default=NULL, help='gene coordinates in reference'),
  make_option(c('-t','--tophits'),     action='store', default=NULL, help='tophit by gene against MEDUSA assembly'),
  make_option(c('-p','--prefix'),      action='store', default=NULL, help='prefix to use in output files'),
  make_option(c('-m','--min'), type="numeric", action='store', default=NULL, help='minimum length to consider a gene')
)

opt <- parse_args(OptionParser(option_list = option_list))
ref = opt$reference
spp = opt$species
genecoords = opt$genecoords
tophits = opt$tophits
prefix = opt$prefix
minlen = opt$min

## read and merge datasets ==============================================================
ref.coord <- read.table(genecoords)
hit.coord <- read.table(tophits, header=T, fill=T)
# filter tophits based on minimum gene length
hit.coord = subset(hit.coord, genelen >= minlen)

## update coordinates in both datasets ==================================================
# merge the two datasets
names(ref.coord) = c("gene", "chrom", "sta", "end", "strand")
mrg <- merge(ref.coord, hit.coord, by="gene", sort=F)
# big chromosome names
if (ref == "hmelv25") { mrg$refchr = substr(mrg$chrom, 6, 7) }
if (ref == "heradem") { mrg$refchr = substr(mrg$chrom, 7, 8) }
if (spp != "hhimfat") { mrg$mapchr1 = substr(mrg$scaffold, 5, 6) }
if (spp == "hhimfat") { mrg$mapchr1 = substr(mrg$scaffold, 8, 9) }
if (spp != "hhimfat") { mrg$mapchr2 = substr(mrg$chr_2ndbest, 5, 6) }
if (spp == "hhimfat") { mrg$mapchr2 = substr(mrg$chr_2ndbest, 8, 9) }

## ==================================================
mrgflt = subset(mrg, refchr == mapchr1)
mrgflt$RefRelPos = 0
mrgflt$QueRelPos = 0
# update reference relative coordinates
df1 = data.frame()
for (bigchr in unique(mrgflt$refchr)) {
  s1 = subset(mrgflt, refchr == bigchr)
  lastPos = 0
  for (scaff in unique(s1$chrom)) {
    s2 = subset(s1, chrom == scaff)
    s2 = s2[order(s2$sta.x),]
    s2$RefRelPos = seq(1,nrow(s2),1) + lastPos
    lastPos = lastPos + nrow(s2)
    df1 = rbind(df1, s2)
  }
}
# update query relative coordinates
df2 = data.frame()
for (bigchr in unique(df1$mapchr1)) {
  s1 = subset(df1, mapchr1 == bigchr)
  lastPos = 0
  for (scaff in unique(s1$chrom)) {
    s2 = subset(s1, chrom == scaff)
    s2 = s2[order(s2$sta.y),]
    s2$QueRelPos = seq(1,nrow(s2),1) + lastPos
    lastPos = lastPos + nrow(s2)
    df2 = rbind(df2, s2)
  }
}
## test correlation of coordinates ==================================================
df.cor = data.frame()
for (ch in unique(df2$refchr)) {
  s1 = subset(df2, refchr == ch)
  cor = cor.test(x=s1$RefRelPos, y=s1$QueRelPos)
  df = data.frame(chrom = ch, rho = round(cor$estimate,2), cor$p.value)
  df.cor = rbind(df.cor, df)
}
write.table(file=paste0(prefix,".coordCorr.PerChromosome.txt"), df.cor, quote=F, row.names=F)

# require(ggplot2)
# ggplot(df2) + geom_point(aes(x=RefRelPos, y=QueRelPos, colour=refchr), alpha=0.2) + facet_wrap(~refchr, ncol=7, scales = "free")
