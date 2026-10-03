nextflow.enable.types = true

process METHURATOR_GTESTIMATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/methurator:2.2.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/methurator:2.2.0--pyhdfd78af_0'}"

    input:
    record(
        meta: Record,
        bam: Path,
        bai: Path,
        fasta: Path,
        fai: Path
    )

    output:
    record(
        meta               : meta,
        methurator_summary : file("${prefix}.yml")
    )

    topic:
    tuple(task.process, 'methurator', eval("methurator --version | sed 's/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    methurator gt-estimator \\
        ${bam} \\
        --fasta ${fasta} \\
        -@ ${task.cpus} \\
        --outdir . \\
        ${args}

    mv methurator_summary.yml ${prefix}.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yml

    """
}
