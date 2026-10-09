

process COMBINED_METHOBJ {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' }"
	
	input:
    path(methylkit_files)
	val(methylkit_meta)

	output:
	path "methobj.rds"       , emit: methobj
	path "methobj_lowcov.rds", emit: methobj_lowcov
	path "versions.yml"      , emit: versions

	script:
	"""

	make_combined_methobj.R \\
		--input_files '$methylkit_files' \\
		--meta '$methylkit_meta'
	

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
	    bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))")
	END_VERSIONS
	"""


}