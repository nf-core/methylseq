

process TWIST_ANNOTATION {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' }"
	
	input:
    path(methobj)
	path(twist_bed_file)

	output:
	path "cpg_obj_twist_out.rds", emit: cpg_obj_twist_out
	path "methobj_lowcov_twist_annotation.rds" , emit: methobj_lowcov_twist_annotation_obj
	path "versions.yml"      , emit: versions

	script:
	"""

	twist_annotation.R \\
		--methobj $methobj \\
		--twist_bed_file $twist_bed_file \\
		--cpg_obj_twist_out cpg_obj_twist_out.rds \\
		--methobj_twist_out methobj_lowcov_twist_annotation.rds

	

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