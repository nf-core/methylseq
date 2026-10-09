// Merge tiled beta values across samples
process MERGE_TILED_STATS {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0'

    input:
    path tiled_methobj

    output:
    path "merge_tiled_stats.csv", emit: merge_tiled_stats
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_tiled_stats.py \\
        --tiled_methobj ${tiled_methobj} \\
        --output merge_tiled_stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))" 2>/dev/null)
    END_VERSIONS
    """
}
