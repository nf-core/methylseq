// Convert methyl_qc outputs into MultiQC custom content
process METHYLQC_MULTIQC {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1'

    input:
    path control_stats
    path read_metrics
    path meth_distribution
    path depth_distribution
    path mbias
    path chr_percent
    path twist_annotation
    path repeat_annotation
    path tiled_beta

    output:
    path "*_mqc.json"  , emit: mqc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    methylqc_multiqc.py \\
        --control_stats ${control_stats} \\
        --read_metrics ${read_metrics} \\
        --meth_distribution ${meth_distribution} \\
        --depth_distribution ${depth_distribution} \\
        --mbias ${mbias} \\
        --chr_percent ${chr_percent} \\
        --twist_annotation ${twist_annotation} \\
        --repeat_annotation ${repeat_annotation} \\
        --tiled_beta ${tiled_beta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methylqc_multiqc: \$(methylqc_multiqc.py --version 2>/dev/null)
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null)
    END_VERSIONS
    """
}
