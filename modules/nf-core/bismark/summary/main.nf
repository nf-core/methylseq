nextflow.preview.types = true

process BISMARK_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38e61d14ccaed82f60c967132963eb467d0fa4bccb7a21404c49b4f377735f03/data' :
        'community.wave.seqera.io/library/bismark:0.25.1--1f50935de5d79c47' }"

    input:
    record(
        bam: Set<String>,
        align_report: Set<Path>,
        dedup_report: Set<Path>,
        methylation_report: Set<Path>,
        methylation_mbias: Set<Path>
    )

    output:
    record(
        html: file("*report.html"),
        txt: file("*report.txt")
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    bismark2summary ${bam.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch bismark_summary_report.txt
    touch bismark_summary_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
