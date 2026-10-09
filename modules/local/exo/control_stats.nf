// Mean CpG methylation of the lambda and pUC19 spike-in controls per sample
process CONTROL_STATS {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0'

    input:
    path methobj

    output:
    path "control_stats.csv", emit: control_stats
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    control_stats.py \\
        --methobj ${methobj} \\
        --stats_out control_stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))" 2>/dev/null)
    END_VERSIONS
    """
}
