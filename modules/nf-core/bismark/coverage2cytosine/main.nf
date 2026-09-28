nextflow.enable.types = true

process BISMARK_COVERAGE2CYTOSINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bddea334e6ccbce005ce540214747acf822b040185d2198220dcfbb4b258c331/data' :
        'community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4' }"

    input:
    record(
        meta: Record,
        methylation_coverage: Path,
        fasta: Path,
        bismark_index: Path
    )

    stage:
    stageAs fasta, 'tmp/*' // This change mounts as directory containing the FASTA file to prevent nested symlinks

    output:
    record(
        id                         : meta.id,
        meta                       : meta,
        coverage2cytosine_coverage : file("*.cov.gz", optional: true),
        coverage2cytosine_report   : file("*report.txt.gz"),
        coverage2cytosine_summary  : file("*cytosine_context_summary.txt")
    )

    topic:
    tuple(task.process, 'bismark', eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage2cytosine \\
        ${methylation_coverage} \\
        --genome ${bismark_index} \\
        --output ${prefix} \\
        --gzip \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.cov.gz
    echo "" | gzip > ${prefix}.report.txt.gz
    touch ${prefix}.cytosine_context_summary.txt
    """
}
