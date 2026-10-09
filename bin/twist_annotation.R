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
	make_option(c("-b", "--twist_bed_file"), type="character", default=NULL, help="TWIST methylome bed file", metavar="character"),
	make_option(c("-c", "--cpg_obj_twist_out"), type="character", default=NULL, help="RDS file with TWIST methylome CpG annotation", metavar="character"),
	make_option(c("-o", "--methobj_twist_out"), type="character", default=NULL, help="Methobj RDS file with TWIST methylome annotation", metavar="character")
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

if (is.null(opt$twist_bed_file)){
    print_help(opt_parser)
    stop("TWIST bed file needs to be provided!")
} else {
    twist_bed_file = opt$twist_bed_file
}

if (is.null(opt$cpg_obj_twist_out)){
    print_help(opt_parser)
    stop("Output name for twist RDS file with TWIST methylome annotation needs to be provided!")
} else {
    cpg_obj_twist_out = opt$cpg_obj_twist_out
}

if (is.null(opt$methobj_twist_out)){
    print_help(opt_parser)
    stop("Output name for Methobj RDS file with TWIST methylome annotation needs to be provided!")
} else {
    methobj_twist_out = opt$methobj_twist_out
}

twist_bed_rds <- function(bed_file_twist){
	methylome_annot_bed = readGeneric(bed_file_twist, meta.cols=4)
	
	# We capture all "generic" annotations in the file with regex filtering.
	all_annot = unique(unlist(lapply(elementMetadata(methylome_annot_bed)[[1]],strsplit,',')))
	annot_filt = all_annot[!grepl('^cg', all_annot)]
	annot_filt = annot_filt[!grepl('^ch', annot_filt)]
	annot_filt = annot_filt[!grepl('^rs[0-9]+', annot_filt)]

	# Create a GRangesList object whose elements are the granges annotated according to each element in annot_filt.
	cpg_obj_twist = c()
	for (annotation_case in annot_filt){
	    this_subset = methylome_annot_bed[grepl(annotation_case, elementMetadata(methylome_annot_bed)[[1]])]
	    cpg_obj_twist = c(cpg_obj_twist, this_subset)
	}
	cpg_obj_twist = GRangesList(cpg_obj_twist)
	names(cpg_obj_twist@partitioning) <- annot_filt
	return(cpg_obj_twist)
}


# Create TWIST methylome object
cpg_obj_twist <- twist_bed_rds(twist_bed_file)
saveRDS(cpg_obj_twist, file=cpg_obj_twist_out)

feat_annot_reps = c()
for (sample in methobj){
	methobj_gr = as(sample,'GRanges')
	feat_annot = c()
	for (annotation_case_index in 1:length(cpg_obj_twist)){
		feat_annot = c(feat_annot, annotateWithFeature(methobj_gr, cpg_obj_twist[annotation_case_index][[1]]))
	}
	# Note: Is this overwriting feat_annot_reps each time in the loop??
	feat_annot_reps = c(feat_annot_reps, feat_annot) 
}
saveRDS(feat_annot_reps, file=methobj_twist_out)
