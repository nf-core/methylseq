/**
This module describes the Rastair process 'mbias' to
run Rastair and assessing C->T conversion
as a readout for methylation
in a per read position basis
*/

process RASTAIR_MBIAS {
    label 'process_medium'
    container "docker.io/sbludwig/rastair:version-0.8.2"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.rastair_mbias.txt"),   emit: txt
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rastair mbias \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        ${bam} > ${prefix}.rastair_mbias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
