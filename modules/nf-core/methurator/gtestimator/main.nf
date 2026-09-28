nextflow.enable.types = true

process METHURATOR_GTESTIMATOR {
    tag id
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/methurator:2.2.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/methurator:2.2.0--pyhdfd78af_0'}"

    input:
    record(
        id: String,
        bam: Path,
        bai: Path,
        fasta: Path,
        fai: Path,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id                 : id,
        methurator_summary : file("${prefix}.yml")
    )

    topic:
    tuple(task.process, 'methurator', eval("methurator --version | sed 's/.* //'")) >> 'versions'

    script:
    args = args ?: ''
    prefix = prefix ?: "${id}"
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
    prefix = prefix ?: "${id}"
    """
    touch ${prefix}.yml

    """
}
