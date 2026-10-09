// Percent of covered CpGs in each Twist category and in repeats, per sample
process MERGE_ANNOTATION_STATS {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0'

    input:
    path methobj
    path methobj_twist_annotation
    path cpg_obj
    path methobj_repeat_annotation

    output:
    path "twist_annotation.csv" , emit: twist_annotation
    path "repeat_annotation.csv", emit: repeat_annotation
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_annotation_stats.py \\
        --methobj ${methobj} \\
        --methobj_twist_annotation ${methobj_twist_annotation} \\
        --cpg_obj ${cpg_obj} \\
        --methobj_repeat_annotation ${methobj_repeat_annotation} \\
        --twist_annotation_csv twist_annotation.csv \\
        --repeat_annotation_csv repeat_annotation.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))" 2>/dev/null)
    END_VERSIONS
    """
}
