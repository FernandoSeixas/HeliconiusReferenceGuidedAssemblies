## load libraries ----
require(optparse)
library(foreach)
library(doMC)
require(dplyr)

## parameters ----
option_list <- list(
	make_option(c('-i','--input'),   action='store', default=NULL, type='character', help='input file'),
	make_option(c('-o','--output'),  action='store', default=NULL, type='character', help='output file'),
	make_option(c('-b','--buffer'),  action='store', default=NULL, type='numeric', help='max. distance btw hits to merge'),
	make_option(c('-c','--cores'),   action='store', default=NULL, type='numeric', help='number of cores')
)

opt <- parse_args(OptionParser(option_list = option_list))
file = opt$input
output=opt$output
buffer = opt$buffer # same gene hits closer then <buffer> are merged and their attributes updated 
cores=opt$cores

## read data and transform ----------
inp <- read.table(file)
inp <- inp[,c(1,5,8,11,12,13,15,20)]
names(inp) = c("gene", "genelen","scaffold","sta","end","identity","maplen","strand")
inp$hits = 1

## merge consecutive windows of gene mapping to multiple regions ------------------------------
print("merging genes")
# list genes in input
genes <- unique(inp$gene)
print(length(genes))
# apply merge function in parallel to genes
source("/n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/mergemaps.R")
registerDoMC(cores=cores)
data <- foreach(n = 1:length(genes), .combine=rbind) %dopar% mergemaps(n)

check1 = length(unique(data$gene))
print(check1)
if (check1 != length(genes)) {stop("not all genes present")}

#  calculate the proportion of the gene that is covered by hits (caution: this is not exact as some hits may overlap)
data$propmap = data$maplen / data$genelen * 100
data$propmap = ifelse(data$propmap > 100, 100, data$propmap) # correct if proportion greater than 100% (obvious artifact of having multiple hits that partially overlap)

# score hits: score = propmap (sum across hits) * percentidentity (average across hits)
data$score = data$propmap * (data$identity/data$hits) / 10000
write.table(file=paste0(output,".mergedmappings.txt"), data, row.names=F, quote=F, sep="\t")


## get top hits --------------------
print("getting tophits")
tophits <- data %>% group_by(gene) %>% top_n(1, score) %>% sample_n(1) %>% as.data.frame
# calculate delta score (difference between 1st and 2nd best hits)
tophits$deltascore <- 999
tophits$chr_2ndbest <- ""
tophits$sta_2ndbest <- -999
tophits$nb.hits <- 0
genes <- unique(tophits$gene)

second_hit <- function(i) {
    if (i %in% seq (100,30000,100)) { return(print(i)) }
    genename = genes[i]
    s1 = subset(data, gene == genename)
    if (nrow(s1) > 1) {
        twohits <- s1 %>% group_by(gene) %>% top_n(2, score) %>% as.data.frame
        tophits$deltascore[i]  <<- abs(twohits$score[1] - twohits$score[2])
	    tophits$chr_2ndbest[i] <<- as.character(twohits$scaffold[2])
	    tophits$sta_2ndbest[i] <<- twohits$sta[2]
    }
}

foreach(n = 1:length(genes)) %do% second_hit(n)
tophits$deltascore = ifelse(tophits$deltascore == 999, 1, tophits$deltascore)

# nb of hits / gene
tophits$nb.hits <- 0
for (rw in 1:nrow(tophits)) {
	genename = tophits$genes[rw]
	tophits$nb.hits = nrow(subset(data, gene == genename))
}

# write to file
write.table(file=paste0(output, ".tophits.txt"), tophits, row.names=F, quote=F, sep="\t")


