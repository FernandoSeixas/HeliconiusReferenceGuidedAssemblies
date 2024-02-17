require(optparse)

# get input parameters
option_list <- list(make_option(c('-g','--gap'), action='store', default=NULL, help='missing sites count file'),
		    make_option(c('-s','--sca'), action='store', default=NULL, help='scaffold size file'),
		    make_option(c('-n','--spp'), action='store', default=NULL, help='species names')
		    )

opt <- parse_args(OptionParser(option_list = option_list))

gap.file <- opt$gap
sca.file <- opt$sca
spp.name <- opt$spp

gap.len <- read.table(gap.file)
sca.len <- read.table(sca.file)
cmb.dat <- cbind(sca.len, gap.len)
names(cmb.dat) <- c("scaff_size", "missN_size")
cmb.dat$perc_miss <- cmb.dat$missN_size / cmb.dat$scaff_size * 100

print(paste(spp.name,sum(cmb.dat$missN_size),sum(cmb.dat$scaff_size),sum(cmb.dat$missN_size)/sum(cmb.dat$scaff_size)))
