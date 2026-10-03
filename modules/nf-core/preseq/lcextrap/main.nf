nextflow.enable.types = true

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_single'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.2.0--hdcf5f25_6':
        'biocontainers/preseq:3.2.0--hdcf5f25_6' }"

    input:
    record(
        meta: Record,
        bam: Path,
        args: String?,
        prefix: String?
    )

    output:
    record(
        meta      : meta,
        lc_extrap : file("*.lc_extrap.txt"),
        lc_log    : file("*.log")
    )

    topic:
    tuple(task.process, 'preseq', eval("preseq 2>&1 | sed -n 's/Version: //p'")) >> 'versions'

    script:
    args = task.ext.args ?: args ?: ''
    args = task.attempt > 1 ? args + ' -defects' : args  // Disable testing for defects
    prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $args \\
        $paired_end \\
        -output ${prefix}.lc_extrap.txt \\
        $bam
    cp .command.err ${prefix}.command.log
    """

    stub:
    prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    """
    touch ${prefix}.lc_extrap.txt
    touch ${prefix}.command.log
    """
}
