process CUSTOM_MULTIQC_PICARDHS {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                '062689759574.dkr.ecr.us-east-1.amazonaws.com/custom_multiqc_picardhs:collect_picardHs_0.0.1' :
                '062689759574.dkr.ecr.us-east-1.amazonaws.com/custom_multiqc_picardhs:collect_picardHs_0.0.1'}"

    input:
    path ("picard/*")
    path sample_sheet
    path picard_bargraph_header
    path picard_linegraph_header
    path picard_table_config

    output:
    path "*_mqc.yml"           , emit: yml
    path "*_mqc.tsv"           , emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    collect_picardHs.py \\
        picard \\
        $sample_sheet \\
        $picard_bargraph_header \\
        $picard_table_config

    cat $picard_linegraph_header <(echo) 'uniformity.yml' > 'picard_uniformity_mqc.yml'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        custom_multiqc_picardHs: \$(echo \$(collect_picardHs.py --version 2>&1) | sed 's/^.*collect_picardHs //; s/Using.*\$//' )
    END_VERSIONS
    """
}
