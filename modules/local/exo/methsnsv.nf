process METHSNSV {
    tag "$meta.id"
    // label 'process_low_memscale'
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                '062689759574.dkr.ecr.us-east-1.amazonaws.com/methsnsv:meth_stats_0.0.1' :
                '062689759574.dkr.ecr.us-east-1.amazonaws.com/methsnsv:meth_stats_0.0.1'}"

    input:
    tuple val(meta), path(bedgraph)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meth_stats.py \\
        $bedgraph \\
        ${meta.id}_snsv.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methsnsv: meth_stats_0.0.1
    END_VERSIONS
    """
}
