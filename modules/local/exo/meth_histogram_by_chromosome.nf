// Per-chromosome CpG counts and percent methylated CpGs
process METH_HISTOGRAM_BY_CHROMOSOME {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1'

    input:
    path bedgraph_files
    val  bedgraph_meta
    path chrom_sizes

    output:
    path "counts_df_long.csv", emit: counts_df_long
    path "chr_percent_df.csv", emit: chr_percent_df
    path "*.png"             , emit: plots
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    meth_histogram_by_chromosome.py \\
        --input_files '${bedgraph_files}' \\
        --meta '${bedgraph_meta}' \\
        --hg38_chrom_sizes ${chrom_sizes} \\
        --counts_df_long counts_df_long.csv \\
        --chr_percent_df chr_percent_df.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null)
    END_VERSIONS
    """
}
