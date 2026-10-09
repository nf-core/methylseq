// Summarise methylation counts in 1 kb tiles
process TILE_METHYL_COUNTS {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0'

    input:
    path methobj_lowcov

    output:
    path "methobj_lowcov_tiled.rds", emit: methobj_lowcov_tiled
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tile_methyl_counts.R \\
        --methobj ${methobj_lowcov} \\
        --window_size 1000 \\
        --step_size 1000 \\
        --cov_bases 0 \\
        --output methobj_lowcov_tiled.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-methylkit: \$(Rscript -e "library(methylKit); cat(as.character(packageVersion('methylKit')))" 2>/dev/null)
    END_VERSIONS
    """
}
