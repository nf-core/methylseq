nextflow.enable.types = true

process RASTAIR_MBIAS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15120636da858ba73a2493281bfa418005f08c0ed09369a837c05f3f9e14a4a6/data' :
        'community.wave.seqera.io/library/rastair:0.8.2--bf70eeab4121509c' }"

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
        meta          : meta,
        rastair_mbias : file("*.rastair_mbias.txt")
    )

    topic:
    tuple(task.process, 'rastair', eval("rastair --version | sed 's/rastair //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rastair mbias \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        ${bam} > ${prefix}.rastair_mbias.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_mbias.txt
    """
}
