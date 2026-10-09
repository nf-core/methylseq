// Per-sample CpG depth and methylation statistics and merged M-bias tables
process BASIC_READ_STATISTICS {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1'

    input:
    path bedgraph_hdf5
    path mbias_files
    val  mbias_meta

    output:
    path "merged_mbias.csv"          , emit: merged_mbias
    path "mbias.id_file.csv"         , emit: mbias_id_file
    path "tab_id_value_counts_df.csv", emit: meth_distribution
    path "rnf_counts_id_df.csv"      , emit: depth_distribution
    path "metrics_table.csv"         , emit: metrics_table
    path "*.png"                     , emit: plots
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    basic_read_statistics.py \\
        --bedGraph_hdf5 ${bedgraph_hdf5} \\
        --mbias_files '${mbias_files}' \\
        --meta '${mbias_meta}' \\
        --id_file mbias.id_file.csv \\
        --merged_mbias merged_mbias.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basic_read_statistics: \$(basic_read_statistics.py --version 2>/dev/null)
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null)
    END_VERSIONS
    """
}
