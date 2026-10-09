

process MERGE_TILED_STATS {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' }"
	
	input:
    path(tiled_methobj)
	val(methylkit_meta)

	output:
	
	path "merge_tiled_stats.csv", emit: merge_tiled_stats
	// path "versions.yml"      , emit: versions

	script:
	"""

	merge_tiled_stats.py \\
		--tiled_methobj $tiled_methobj \\
		--meta '$methylkit_meta' \\
		--output merge_tiled_stats.csv
	

	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))")
    END_VERSIONS
	"""


}