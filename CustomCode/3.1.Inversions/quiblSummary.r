##### load libraries /////////////////////////////////////////////////////
library("quiblR")
library("ggplot2")
library("ape")
library("hash")
library("ggtree")
library("ggpubr")
library("dplyr")



##### read data //////////////////////////////////////////////////////////
win = "10kb_50kb"
flt = ".minSites500"
emr = ".EM50"
set = paste0(win,flt,emr)
dir = "/n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/4.1.Inversions/minimap2/quibl/"
speciesTree <- read_speciesTree(paste0(dir,"speciesTree.nwk"))
quiblOutput <- read_csv_quibl(paste0(dir,"chromosomes.phyml_bionj.",set,".quibl.csv"))
originalTrees <- read.tree(paste0(dir,"chromosomes.phyml_bionj.",paste0(win,flt),".trees"))



##### Get some big-picture results ////////////////////////////////////////
quiblOutput <- mutate(quiblOutput, isDiscordant=as.integer(! apply(quiblOutput, 1, isSpeciesTree, sTree=speciesTree)))
quiblOutput <- mutate(quiblOutput, isSignificant = as.numeric(apply(quiblOutput, 1, testSignificance, threshold=10)))
totalTrees <- sum(quiblOutput$count)/length(unique(quiblOutput$triplet))
quiblOutput <- mutate(quiblOutput,totalIntrogressionFraction=(mixprop2*count*isDiscordant)/totalTrees)
quiblOutput



##### Evaluating model fit - does it pass the smell test ////////////////////////////////////////

# get triplet of interest
tripletsList = unique(quiblOutput$triplet)
TargetTriplet = tripletsList[grepl("burA", tripletsList)]
# TargetTriplet = tripletsList[grepl("dorA", tripletsList)]
TargetTriplet = TargetTriplet[grepl("melA", TargetTriplet)]
TargetTriplet = TargetTriplet[grepl("eraA", TargetTriplet)]

# define outgroup of tripled
outg = "melA"

# per locus statistics 
perLocusStats <- getPerLocusStats(quiblOutput = quiblOutput, trip = TargetTriplet, treeList = originalTrees, overallOut="eueA")
head(perLocusStats)

# calculate expected distrubutions
introTop_ILSOnly <- getILSOnlyDist(0,0.1,subset(quiblOutput, outgroup==outg & triplet==TargetTriplet))
introTop_ILSMix <- getILSMixtureDist(0,0.1,subset(quiblOutput, outgroup==outg & triplet==TargetTriplet))
introTop_nonILSMix <- getNonILSMixtureDist(0,0.1,subset(quiblOutput, outgroup==outg & triplet==TargetTriplet))

# combine expected distributions into data.frame
introComb = data.frame(t2 = introTop_ILSOnly$x, 
                       ILSOnly = introTop_ILSOnly$y,
                       ILSMix = introTop_ILSMix$y,
                       nonILSMix = introTop_nonILSMix$y
)
# plot observed + expected distributions
subData = subset(perLocusStats, out == outg)
# ggplot() +
#   geom_histogram(aes(subData$branchLength), binwidth = 0.001, fill="black", col="gray30", alpha=0.5) +
#   # geom_line(data=introComb, aes(x=t2, y=ILSOnly), col="blue") +
#   geom_line(data=introComb, aes(x=t2, y=nonILSMix+ILSMix), col="red") +
#   geom_line(data=introComb, aes(x=t2, y=ILSMix), col="blue") +
#   geom_line(data=introComb, aes(x=t2, y=nonILSMix), col="cyan")





##### TARGET REGION #####

## get trees from target region
data = read.table(paste0(dir,"chromosomes.phyml_bionj.",paste0(win,flt),".data.tsv"))
names(data) = c("scaffold","start","end","midSite", "sites", "-lnL")
perLocusStatsExt = cbind(perLocusStats, data)
targetRegion = subset(perLocusStatsExt, scaffold == "Hmel213001o" & start >= 16917575 & end <= 17154514)
# targetRegion = subset(perLocusStatsExt, scaffold == "Hmel213001o" & start <= 17154514 & end >= 16917575)

## Introgression probability vs T2 branch length

# calculate introrgession probability as function of t2 branch length
introComb$probI = introComb$nonILSMix/(introComb$nonILSMix+introComb$ILSMix)

coeff = 200
a = ggplot() +
  geom_histogram(aes(subData$branchLength), binwidth = 0.001, fill="black", col="gray30", alpha=0.4) +
  # geom_line(data=introComb, aes(x=t2, y=ILSOnly), lty=1) +
  geom_line(data=introComb, aes(x=t2, y=ILSMix), lty=2) +
  geom_line(data=introComb, aes(x=t2, y=nonILSMix), lty=3) +
  geom_vline(xintercept = mean(targetRegion$branchLength), col="blue") +
  # geom_vline(data=targetRegion, aes(xintercept = branchLength), col="blue") +
  geom_line(data = introComb, aes(x=t2, y=probI*coeff)) +
  xlim(0,0.050) +
  scale_y_continuous(
    # first axis 
    name = "Count",
    # add secondary axis
    sec.axis = sec_axis( trans= ~./coeff, name="Introgression Probability")
  )  +
  xlab("Internal Branch Length") +
  ggtitle("H. burneyi") 


require(gridExtra)
grid.arrange(a, b, nrow=2)
