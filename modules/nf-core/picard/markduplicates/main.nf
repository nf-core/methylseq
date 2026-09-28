nextflow.enable.types = true

process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    record(
        meta: Record,
        bam: Path,
        fasta: Path?,
        fai: Path?
    )

    output:
    record(
        meta           : meta,
        bam            : file("*.bam", optional: true),
        bai            : file("*.bai", optional: true),
        cram           : file("*.cram", optional: true),
        picard_metrics : file("*.metrics.txt")
    )

    topic:
    tuple(task.process, 'picard', eval("picard MarkDuplicates --version 2>&1 | sed -n 's/.*Version://p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${bam.getExtension()}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    if ("${bam}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.${suffix} \\
        ${reference} \\
        --METRICS_FILE ${prefix}.metrics.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${bam.getExtension()}"
    if ("${bam}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.${suffix}.bai
    touch ${prefix}.metrics.txt
    """
}
