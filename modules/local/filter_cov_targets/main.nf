// Keep the Bismark coverage rows (1-based `pos pos`) whose position falls in a target, each row once and unchanged
process FILTER_COV_TARGETS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(cov), path(targets)

    output:
    tuple val(meta), path("*.targeted.cov"), emit: intersect
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${cov.baseName}"
    """
    # Each row becomes [pos-1, pos) followed by the original row, which is written back unchanged
    case "${cov}" in
        *.gz) zcat ${cov} ;;
        *) cat ${cov} ;;
    esac \\
        | awk 'BEGIN{FS=OFS="\t"} NF{print \$1, \$2-1, \$2, \$0}' \\
        | bedtools intersect -a stdin -b ${targets} -wa -u \\
        | cut -f 4- \\
        > ${prefix}.targeted.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${cov.baseName}"
    """
    touch ${prefix}.targeted.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
