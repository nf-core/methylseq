// Annotate low-coverage methylation calls with Twist methylome target categories
process TWIST_ANNOTATION {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0'

    input:
    path methobj
    path twist_bed_file

    output:
    path "cpg_obj_twist_out.rds"              , emit: cpg_obj_twist_out
    path "methobj_lowcov_twist_annotation.rds", emit: methobj_lowcov_twist_annotation_obj
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    twist_annotation.R \\
        --methobj ${methobj} \\
        --twist_bed_file ${twist_bed_file} \\
        --cpg_obj_twist_out cpg_obj_twist_out.rds \\
        --methobj_twist_out methobj_lowcov_twist_annotation.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))" 2>/dev/null)
        bioconductor-genomation: \$(Rscript -e "library(genomation); cat(as.character(packageVersion('genomation')))" 2>/dev/null)
        bioconductor-rtracklayer: \$(Rscript -e "library(rtracklayer); cat(as.character(packageVersion('rtracklayer')))" 2>/dev/null)
        bioconductor-chippeakanno: \$(Rscript -e "library(ChIPpeakAnno); cat(as.character(packageVersion('ChIPpeakAnno')))" 2>/dev/null)
        bioconductor-txdb.hsapiens.ucsc.hg38.knowngene: \$(Rscript -e "library(TxDb.Hsapiens.UCSC.hg38.knownGene); cat(as.character(packageVersion('TxDb.Hsapiens.UCSC.hg38.knownGene')))" 2>/dev/null)
    END_VERSIONS
    """
}
