nextflow.preview.types = true

process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7' :
        'biocontainers/methyldackel:0.6.1--he4a0461_7' }"

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
        id: meta.id,
        meta: meta,
        methyldackel_mbias: file("*.mbias.txt")
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MethylDackel mbias \\
        $args \\
        $fasta \\
        $bam \\
        $prefix \\
        --txt \\
        > ${prefix}.mbias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mbias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
