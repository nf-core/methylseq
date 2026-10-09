process BEDTOOLS_PERBASE_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(intervals)
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: genomecov
    path  "versions.yml"                   , emit: version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_list = args.tokenize()

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        genomecov \\
        -dz \\
        -ibam $intervals \\
        $args \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
