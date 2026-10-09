#!/usr/bin/env Rscript

invisible( lapply(c(
    "optparse",
	"methylKit"
), library, character.only=T))

option_list = list(
    make_option(c("-i", "--input_files"), type="character", default=NULL, help="List of input methylkit files (i.e. 'healthy_1.methylkit healthy_2.methylkit disease_1.methylkit disease_2.methylkit'", metavar="character"),
    make_option(c("-m", "--meta"), type="character", default=NULL, help="Sample ID and group information metadata (i.e. [[id:healthy_1, group:healthy], [id:healthy_2, group:healthy], [id:disease_1, group:disease], [id:disease_2, group:disease]]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate and read input
if (is.null(opt$input_files)){
    print_help(opt_parser)
    stop("Input methylkit files needs to be provided!")
} else {
    input_files_string = opt$input_files
}

if (is.null(opt$meta)){
    print_help(opt_parser)
    stop("meta data string needs to be provided!")
} else {
    meta_string = opt$meta
}

get_file_list <- function(file_string){

	# Removing square brackets from the string
	cleaned_string = gsub("^\\[|\\]$", "", file_string)	

	# Splitting the string into element[s
	file_vector = unlist(strsplit(cleaned_string, " \\s*")) # delimit by space
	
	return(as.list(file_vector))
}

char2intVec <- function(charVec){
	# Find unique values and sort them alphabetically
	unique_values = sort(unique(charVec))

	# Assign numbers based on alphabetical rank
	converted_vector = match(charVec, unique_values) - 1

	return(converted_vector)
}


methylkit_files <- get_file_list(input_files_string)

# Extracting id values
id_values <- regmatches(meta_string, gregexpr("(?<=id:)[^,\\]]+", meta_string, perl = TRUE))[[1]]

# Extracting group values
group_values <- regmatches(meta_string, gregexpr("(?<=group:)[^,\\]]+", meta_string, perl = TRUE))[[1]]

# Convert group character names to integers (i.e. c("healthy", "healthy", "disease", "disease") --> c(1,1,0,0))
group_values_int <- char2intVec(group_values)

# Samplesheets without a group column get a single treatment group
if (length(group_values_int) == 0) {
	group_values_int <- rep(0, length(id_values))
}

methobj <- methRead(methylkit_files,
            sample.id = as.list(id_values),
            assembly="hg38",
            treatment=group_values_int,
            context="CpG",
            mincov = 10
)

saveRDS(methobj, file="methobj.rds")

methobj_lowcov <- methRead(methylkit_files,
            sample.id = as.list(id_values),
            assembly="hg38",
            treatment=group_values_int,
            context="CpG",
            mincov = 3
)

saveRDS(methobj_lowcov, file="methobj_lowcov.rds")
