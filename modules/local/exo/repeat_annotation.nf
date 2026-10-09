

process REPEAT_ANNOTATION {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' }"
	
	input:
    path(methobj)
	path(repeat_annot)

	output:
	path "methobj_lowcov_repeat_annotation.rds" , emit: methobj_lowcov_repeat_annotation_obj
	path "versions.yml"      , emit: versions

	script:
	"""

	repeat_annotation.R \\
		--methobj $methobj \\
		--repeat_annot $repeat_annot \\
		--methobj_repeat_out methobj_lowcov_repeat_annotation.rds

	

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
	    bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))")
		bioconductor-genomation: \$(Rscript -e "library(genomation); cat(as.character(packageVersion('genomation')))")
		bioconductor-ChIPpeakAnno: \$(Rscript -e "library(ChIPpeakAnno); cat(as.character(packageVersion('ChIPpeakAnno')))")
		bioconductor-TxDb.Hsapiens.UCSC.hg38.knownGene: \$(Rscript -e "library(TxDb.Hsapiens.UCSC.hg38.knownGene); cat(as.character(packageVersion('TxDb.Hsapiens.UCSC.hg38.knownGene')))")
		bioconductor-rtracklayer: \$(Rscript -e "library(rtracklayer); cat(as.character(packageVersion('rtracklayer')))")
	END_VERSIONS
	"""


}