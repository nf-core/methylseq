nextflow.enable.types = true

process BISMARK_DEDUPLICATE {
    tag id
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bddea334e6ccbce005ce540214747acf822b040185d2198220dcfbb4b258c331/data' :
        'community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4' }"

    input:
    record(
        id: String,
        single_end: Boolean,
        bam: Path,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id           : id,
        bam          : file("*.deduplicated.bam"),
        dedup_report : file("*.deduplication_report.txt")
    )

    topic:
    tuple(task.process, 'bismark', eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'")) >> 'versions'

    script:
    args = args ?: ''
    def seqtype = single_end ? '-s' : '-p'
    """
    deduplicate_bismark \\
        ${args} \\
        ${seqtype} \\
        --bam ${bam}
    """

    stub:
    prefix = prefix ?: "${id}"
    """
    touch ${prefix}.deduplicated.bam
    touch ${prefix}.deduplication_report.txt
    """
}
