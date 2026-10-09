#!/usr/bin/env Rscript

invisible( lapply(c(
    "optparse",
	"methylKit"
), library, character.only=T))


option_list = list(
    make_option(c("-m", "--methobj"), type="character", default=NULL, help="Methylkit object RDS file", metavar="character"),
    make_option(c("-w", "--window_size"), type="integer", default=1000, help="Window size for tiling methylkit counts", metavar="character"),
	make_option(c("-s", "--step_size"), type="integer", default=1000, help="Step size for tiling methylkit counts", metavar="character"),
	make_option(c("-c", "--cov_bases"), type="integer", default=0, help="Cov_bases for tiling methylkit counts", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, help="Output tiled Methylkit object RDS file", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate and read input
if (is.null(opt$methobj)){
    print_help(opt_parser)
    stop("Input methylkit object needs to be provided!")
} else {
    methobj = readRDS(opt$methobj)
}
if (is.null(opt$output)){
    print_help(opt_parser)
    stop("Input methylkit object needs to be provided!")
} else {
    outputRDS = opt$output
}

if (is.null(opt$window_size)){
    print_help(opt_parser)
    stop("Window size needs to be provided!")
} else {
    window_size = opt$window_size
}

if (is.null(opt$step_size)){
    print_help(opt_parser)
    stop("Window size needs to be provided!")
} else {
    step_size = opt$step_size
}

if (is.null(opt$cov_bases)){
    print_help(opt_parser)
    stop("Window size needs to be provided!")
} else {
    cov_bases = opt$cov_bases
}


methobj_tiles <- tileMethylCounts(
	methobj,
	win_size = window_size, 
	step_size = step_size, 
	cov_bases = cov_bases
)

saveRDS(methobj_tiles,outputRDS)