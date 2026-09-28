nextflow.enable.types = true

process SAMTOOLS_STATS {
    tag id
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    record(
        id: String,
        bam: Path,
        bai: Path,
        fasta: Path?,
        fai: Path?,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id             : id,
        samtools_stats : file("*.stats")
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools version | sed "1!d;s/.* //"')) >> 'versions'

    script:
    args = args ?: ''
    prefix = prefix ?: "${id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${bam} \\
        > ${prefix}.stats
    """

    stub:
    prefix = prefix ?: "${id}"
    """
    touch ${prefix}.stats
    """
}
