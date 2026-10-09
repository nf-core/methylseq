#!/usr/bin/env Rscript

invisible( lapply(c(
    "optparse",
	"methylKit",
	"genomation",
	"ChIPpeakAnno",
	"TxDb.Hsapiens.UCSC.hg38.knownGene",
	"rtracklayer"
), library, character.only=T))






option_list = list(
    make_option(c("-i", "--methobj"), type="character", default=NULL, help="Methobj RDS file", metavar="character"),
	make_option(c("-b", "--repeat_annot"), type="character", default=NULL, help="Repeat annotation file", metavar="character"),
	make_option(c("-o", "--methobj_repeat_out"), type="character", default=NULL, help="Methobj RDS file with TWIST methylome annotation", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate and read input
if (is.null(opt$methobj)){
    print_help(opt_parser)
    stop("Input methobjs files needs to be provided!")
} else {
    methobj = readRDS(opt$methobj)
}

if (is.null(opt$repeat_annot)){
    print_help(opt_parser)
    stop("Repeat annotation file needs to be provided!")
} else {
    repeat_annot = import.bb(opt$repeat_annot)
}

if (is.null(opt$methobj_repeat_out)){
    print_help(opt_parser)
    stop("Output name for Methobj RDS file with TWIST methylome annotation needs to be provided!")
} else {
    methobj_repeat_out = opt$methobj_repeat_out
}


feat_annot_reps = c()
for (sample in methobj){
    methobj_gr = as(sample,'GRanges')
    feat_annot = annotateWithFeature(methobj_gr, repeat_annot)
    feat_annot_reps = c(feat_annot_reps, feat_annot)
}
saveRDS(feat_annot_reps, file=methobj_repeat_out)

