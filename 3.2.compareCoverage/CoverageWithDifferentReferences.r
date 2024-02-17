## load modules
require(stringr)
require(ggplot2)
require(gridExtra)


######################### hec/ele/par Coverage to hmelv25 reference #########################
Species = "HEC|ELE|PAR"
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.2.compareMap/1.2.rmdup/hmelv25/"
lfiles = list.files(dir, pattern = "cov.w25s25.txt")
lfiles = lfiles[grepl(Species, lfiles)]

trio2hmelv25.Coverage = data.frame()
for (file in lfiles) {
  inp = read.table(paste0(dir,file))
  names(inp) = c("scaffold","sta","end","sites","Coverage","Species")
  inp$relCoverage = inp$Coverage / median(inp$Coverage)
  trio2hmelv25.Coverage = rbind(trio2hmelv25.Coverage, inp)
}


######################### hec/ele/par Coverage to MEDUSA reference #########################
Species = "HEC|ELE|PAR"
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.2.compareMap/1.2.rmdup/DISCmel/"
lfiles = list.files(dir, pattern = "cov.w25s25.txt")
lfiles = lfiles[grepl(Species, lfiles)]

trio2medusa.Coverage = data.frame()
for (file in lfiles) {
  inp = read.table(paste0(dir,file))
  names(inp) = c("scaffold","sta","end","sites","Coverage","Species")
  inp$relCoverage = inp$Coverage / median(inp$Coverage)
  trio2medusa.Coverage = rbind(trio2medusa.Coverage, inp)
}


######################### melpomene Coverage to MEDUSA reference #########################
Species = "hec|ele|par"
dir = "/n/mallet_lab/fseixas/1.projects/1.pseudo_references/4.3.GenomeExpansion/mel2MEDUSA/1.2.rmdup/"
lfiles = list.files(dir, pattern = "cov.w25s25.txt")
lfiles = lfiles[grepl(Species, lfiles)]

mel2medusa.Coverage = data.frame()
for (file in lfiles) {
  inp = read.table(paste0(dir,file))
  names(inp) = c("scaffold","sta","end","sites","Coverage","Species")
  inp$relCoverage = inp$Coverage / median(inp$Coverage)
  mel2medusa.Coverage = rbind(mel2medusa.Coverage, inp)
}



##### Transform Datasets as necessary 
trio2hmelv25.Coverage$set = "trio2hmel25"
trio2medusa.Coverage$set = "trio2medusa"
mel2medusa.Coverage$set = "hmel2medusa"

# rename species
trio2hmelv25.Coverage$Species = tolower(trio2hmelv25.Coverage$Species)
trio2medusa.Coverage$Species = tolower(trio2medusa.Coverage$Species)
mel2medusa.Coverage$Species = tolower(mel2medusa.Coverage$Species)

trio2hmelv25.Coverage$Species = str_remove(trio2hmelv25.Coverage$Species, "el_")
trio2medusa.Coverage$Species = str_remove(trio2medusa.Coverage$Species, "el_")
mel2medusa.Coverage$Species = str_remove(mel2medusa.Coverage$Species, "el_")

trio2hmelv25.Coverage$Species = factor(trio2hmelv25.Coverage$Species, levels = c("hhec","hele","hpar"))
trio2medusa.Coverage$Species = factor(trio2medusa.Coverage$Species, levels = c("hhec","hele","hpar"))

# add reference for the mel2medusa dataset
mel2medusa.Coverage$Reference = substr(mel2medusa.Coverage$scaffold, 1,4)
mel2medusa.Coverage$Reference = factor(mel2medusa.Coverage$Reference, levels = c("hhec","hele","hpar"))

# chromosomes
bigchroms = aggregate(trio2medusa.Coverage$end, list(trio2medusa.Coverage$scaffold), max)
bigchroms = subset(bigchroms, x > 8000000)[,1]
hhecBigChroms = bigchroms[grepl("hhec", bigchroms)]
heleBigChroms = bigchroms[grepl("hele", bigchroms)]
hparBigChroms = bigchroms[grepl("hpar", bigchroms)]
#/
hmelBigChroms = aggregate(trio2hmelv25.Coverage$end, list(trio2hmelv25.Coverage$scaffold), max)
hmelBigChroms = subset(hmelBigChroms, x > 8000000)[,1]

# subset to chromosome
ch = 8
melChr = as.character(hmelBigChroms[ch])
hecChr = as.character(hhecBigChroms[ch])
eleChr = as.character(heleBigChroms[ch])
parChr = as.character(hparBigChroms[ch])

tri2mel = subset(trio2hmelv25.Coverage, scaffold == as.character(hmelBigChroms[ch]))
tri2med = subset(trio2medusa.Coverage, scaffold %in% c(hecChr,eleChr,parChr) )
mel2med = subset(mel2medusa.Coverage, scaffold %in% c(hecChr,eleChr,parChr) )

# define repeat region limits
# if (melChr == "Hmel202001o") { lim1 = 4075000; lim2 = 4125000 }
# if (melChr == "Hmel208001o") { lim1 = 3300000; lim2 = 3475000 }
# if (melChr == "Hmel204001o") { lim1 = 5650000; lim2 = 5875000 }
# if (melChr == "Hmel209001o") { lim1 = 5125000; lim2 = 5450000 }


# plots
p1 = ggplot(tri2mel) +
  geom_line(aes(x=(sta+end)/2000000, y=relCoverage, col=Species)) +
  geom_hline(yintercept = 1, lty=2) + 
  scale_x_continuous(breaks = seq(0,30,2)) +
  facet_wrap(~Species, ncol=3, scales="free_x") +
  scale_color_manual(values = c("#1b9e77","#7570b3","#d95f02")) +
  xlab("") + ylab("Relative Coverage") + ggtitle("A") +
  theme(legend.position = "none")
# p1 = p1 + annotate("rect", xmin=lim1/1000000, xmax=lim2/1000000, ymin=0, ymax=max(tri2mel$relCoverage), alpha=1.0, fill="grey")
p2 = ggplot(tri2med) +
  geom_line(aes(x=(sta+end)/2000000, y=relCoverage, col=Species)) +
  geom_hline(yintercept = 1, lty=2) +
  scale_x_continuous(breaks = seq(0,30,2)) +
  facet_wrap(~Species, ncol=3, scales="free_x") +
  scale_color_manual(values = c("#1b9e77","#7570b3","#d95f02")) +
  xlab("") + ylab("Relative Coverage") + ggtitle("B") +
  theme(legend.position = "none")
p3 = ggplot(mel2med) +
  geom_line(aes(x=(sta+end)/2000000, y=relCoverage, col=Reference)) +
  geom_hline(yintercept = 1, lty=2) +
  scale_x_continuous(breaks = seq(0,30,2)) +
  facet_wrap(~Reference, ncol=3, scales="free_x") +
  scale_color_manual(values = c("#1b9e77","#7570b3","#d95f02")) +
  xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle("C") +
  theme(legend.position = "none")


g = grid.arrange(p1, p2, p3, nrow=3)
# ggsave(paste0(dir,"chr",ch,".RelativeCoverages.png"), g)




##### subset repeat regions (in the MEDUSA assemblie) /////////////////////////
head(mel2med)

hhecRR = subset(mel2med, Reference == "hhec" & relCoverage < 0.25)
heleRR = subset(mel2med, Reference == "hele" & relCoverage < 0.25)
hparRR = subset(mel2med, Reference == "hpar" & relCoverage < 0.25)

source("~/code/heliconius_seixas/0.general/Rfunctions.R")
mergeRegions(hhecRR[,1:3], 100000)
mergeRegions(heleRR[,1:3], 100000)
mergeRegions(hparRR[,1:3], 100000)

hhec040001-3900000-5225000
hele040001-4175000-5150000
hpar040001-4125000-5250000

hhec020001 4525000 4875000
hele020001 4900000 5275000
hpar020001 4800000 5150000
# ggplot(mel2med) +
#   geom_line(aes(x=(sta+end)/2000000, y=relCoverage, col=Reference)) +
#   geom_hline(yintercept = 1, lty=2) +
#   facet_wrap(~Reference, nrow=3, scales="free_x") +
#   scale_color_manual(values = c("#1b9e77","#7570b3","#d95f02")) +
#   scale_x_continuous(breaks = seq(0,30,0.5)) +
#   xlab("Chromosome Position (Mb)") + ylab("Relative Coverage") + ggtitle("C") +
#   theme(legend.position = "none") + ylim(0,2.0)
# 
# source("/n/home12/fseixas/code/heliconius_seixas/0.general/Rfunctions.R")
# mergeRegions(subset(mel2med, relCoverage < 0.5 & Reference == "hhec"), 100000)
# mergeRegions(subset(mel2med, relCoverage < 0.5 & Reference == "hele"), 100000)
# mergeRegions(subset(mel2med, relCoverage < 0.5 & Reference == "hpar"), 100000)
# 
# 
# 
# 
# (8925000-5575000)/1000000
# (9625000-5700000)/1000000
# (9625000-5500000)/1000000
# 
# 
