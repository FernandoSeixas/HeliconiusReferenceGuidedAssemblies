# load libraries ----
require(optparse)
require(dplyr)
require(reshape2)

## parameters ----
option_list <- list(
  make_option(c('-r','--reference'),  action='store', default=NULL, help='reference'),
  make_option(c('-s','--species'),    action='store', default=NULL, help='species'),
  make_option(c('-g','--genecoords'), action='store', default=NULL, help='gene coordinates in reference'),
  make_option(c('-t','--tophits'),    action='store', default=NULL, help='tophit by gene against MEDUSA assembly'),
  make_option(c('-p','--prefix'),     action='store', default=NULL, help='prefix for output files'),
  make_option(c('-m','--min'), type="numeric", action='store', default=NULL, help='minimum length to consider a gene')
)

opt <- parse_args(OptionParser(option_list = option_list))
ref = opt$reference
spp = opt$species
genecoords = opt$genecoords
tophits = opt$tophits
prefix = opt$prefix
minlen = opt$min


## ---------------------------------------------------------------------
## read and merge datasets ---------------------------------------------
ref.coord <- read.table(genecoords)
hit.coord <- read.table(tophits, header=T, fill=T)
names(ref.coord) = c("gene", "chrom", "sta", "end", "strand")
# subset genes based on minimum length
hit.coord = subset(hit.coord, genelen >= minlen)
# merge the two datasets
mrg <- merge(ref.coord, hit.coord, by="gene", sort=F)
mrg$chr_2ndbest = ifelse(mrg$chr_2ndbest == "-999", "NA", as.character(mrg$chr_2ndbest))

if (ref == "hmelv25") { mrg$refchr = substr(mrg$chrom, 6, 7) }
if (ref == "heradem") { mrg$refchr = substr(mrg$chrom, 7, 8) }
if (spp != "hhimfat") { mrg$mapchr1 = substr(mrg$scaffold, 5, 6) }
if (spp == "hhimfat") { mrg$mapchr1 = substr(mrg$scaffold, 8, 9) }
if (spp != "hhimfat") { mrg$mapchr2 = substr(mrg$chr_2ndbest, 5, 6) }
if (spp == "hhimfat") { mrg$mapchr2 = substr(mrg$chr_2ndbest, 8, 9) }



## get metric of correct scaffold-chromosoe hits -----------------------------------------------------
## check big chromosomes [for which correlations can reasonably be assessed]
big.chroms = paste0(sprintf("%02d", seq(1,21)))
mrg.flt = mrg[mrg$refchr %in% big.chroms,] # consider only hits to "big" chromosomes
mrg.flt = mrg.flt[,c(20,21,22)]

# count correct and incorrect hits (and among the latter how many times the 2nd best is correct)
nb1G = 0
nb2G = 0
allB = 0
hits.df = data.frame()
prop.df = data.frame()
ch = "01"
for (ch in unique(mrg.flt$refchr)) {
  print(ch)
  s1 = subset(mrg.flt, refchr == ch)
  sgood = s1[s1$mapchr1 == ch,]
  sbads = s1[s1$mapchr1 != ch,]
  sbads_good = sbads[sbads$mapchr2 == ch,]

  df1 = data.frame(chrom = ch, nb1_correct = nrow(sgood)/nrow(s1)*100, nb2_correct = nrow(sbads_good)/nrow(s1)*100, nb1_incorrect = (nrow(sbads)-nrow(sbads_good))/nrow(s1)*100, total=nrow(s1)) 
  df2 = data.frame(chrom = ch, nb1_correct = nrow(sgood), nb2_correct = nrow(sbads_good), nb1_incorrect = (nrow(sbads)-nrow(sbads_good)), total=nrow(s1))
  nb1G = nb1G + nrow(sgood); nb2G = nb2G + nrow(sbads_good); allB = allB + (nrow(sbads)-nrow(sbads_good))
  hits.df = rbind(hits.df, df2)
  prop.df = rbind(prop.df, df1)
}
write.table(file=paste0(prefix, ".classHits.txt"), hits.df, quote=F, row.names=F, sep="\t")


## rename hit types and assign position for graph
prop.df.long = melt(prop.df, id.vars = c("chrom"))
names(prop.df.long) = c("chrom", "hittype", "prop")

prop.df.long$hittype = ifelse(
  prop.df.long$hittype == "nb1_correct", "nb1",
  ifelse(prop.df.long$hittype == "nb1_incorrect", "nb3", "nb2")
)

prop.df.long$lab = ifelse(
  prop.df.long$hittype == "nb1", paste0(round(prop.df.long$prop),"%"), 
  ifelse(prop.df.long$hittype == "nb3", paste0(round(prop.df.long$prop),"%"), "")
)

prop.df.long$pos = ifelse(
  prop.df.long$hittype == "nb1", 90, 
  ifelse(prop.df.long$hittype == "nb3", 3, 50)
)

# plot summary of hits
require(ggplot2)
require(ggthemes)
png(paste0(prefix, ".propCorrectHits.perChrom.png"))
p <- ggplot(prop.df.long, aes(y=prop, x=chrom, fill=hittype, label=lab)) +
  geom_bar(stat="identity") +
  geom_text(colour="white", size=4, aes(y=pos)) +
  scale_y_continuous(breaks = seq(0,100,10)) +
  coord_flip() + theme_economist() + scale_fill_economist() +
  xlab("Hmel2.5 Chromosome") + ylab("Gene hits (%)")
plot(p)
dev.off()

nb1_prop = subset(prop.df.long, hittype == "nb1")
write.table(file=paste0(prefix, ".nb1_propHits.txt"), nb1_prop, quote=F, row.names=F, sep="\t")


