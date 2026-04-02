nextflow.preview.types = true

process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    record(
        meta: Record,
        reads: Path,
        fasta: Path,
        fai: Path
    )

    output:
    record(
        id             : meta.id,
        meta           : meta,
        bam            : file("*.bam", optional: true),
        bai            : file("*.bai", optional: true),
        cram           : file("*.cram", optional: true),
        picard_metrics : file("*.metrics.txt"),
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $reads \\
        --OUTPUT ${prefix}.${suffix} \\
        $reference \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
