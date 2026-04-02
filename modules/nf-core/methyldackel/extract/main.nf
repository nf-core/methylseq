nextflow.preview.types = true

process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

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
        id                    : meta.id,
        meta                  : meta,
        methydackel_bedgraph  : file("*.bedGraph", optional: true),
        methydackel_methylkit : file("*.methylKit", optional: true)
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    MethylDackel extract \\
        $args \\
        $fasta \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def out_extension = args.contains('--methylKit') ? 'methylKit' : 'bedGraph'
    """
    touch ${bam.baseName}_CpG.${out_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
