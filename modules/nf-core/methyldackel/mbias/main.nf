nextflow.enable.types = true

process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7' :
        'quay.io/biocontainers/methyldackel:0.6.1--he4a0461_7' }"

    input:
    record(
        meta: Record,
        bam: Path,
        bai: Path,
        fasta: Path,
        fai: Path,
        args: String?,
        prefix: String?
    )

    output:
    record(
        meta               : meta,
        methyldackel_mbias : file("*.mbias.txt")
    )

    topic:
    tuple(task.process, 'methyldackel', eval("MethylDackel --version 2>&1 | cut -f1 -d' '")) >> 'versions'

    script:
    args = task.ext.args ?: args ?: ''
    prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    """
    MethylDackel mbias \\
        ${args} \\
        ${fasta} \\
        ${bam} \\
        ${prefix} \\
        --txt \\
        > ${prefix}.mbias.txt
    """

    stub:
    prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    """
    touch ${prefix}.mbias.txt
    """
}
