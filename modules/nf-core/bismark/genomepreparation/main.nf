process BISMARK_GENOMEPREPARATION {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::bismark=0.24.0"

    input:
    path fasta, stageAs: "BismarkIndex/*"

    output:
    path "BismarkIndex" , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark_genome_preparation \\
        $args \\
        BismarkIndex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
