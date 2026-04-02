nextflow.preview.types = true

process BISMARK_COVERAGE2CYTOSINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38e61d14ccaed82f60c967132963eb467d0fa4bccb7a21404c49b4f377735f03/data' :
        'community.wave.seqera.io/library/bismark:0.25.1--1f50935de5d79c47' }"

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
    file("versions.yml") >> 'versions'

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.cov.gz
    echo | gzip > ${prefix}.report.txt.gz
    touch ${prefix}.cytosine_context_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
