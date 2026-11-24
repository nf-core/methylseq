/*
 * Summarise MethylDackel bedGraphs, Bismark coverage files or rastair 0.8 calls as MultiQC custom content:
 * General Statistics columns, a per-CpG methylation distribution and spike-in/control methylation.
 */
process METHYLATION_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/759bb68f05a46782e9715e821380f641adfb2ba92d41a0b7bfc98536086e101c/data'
        : 'community.wave.seqera.io/library/gawk_gzip:f7d38a45baa28337'}"

    input:
    tuple val(meta), path(calls, stageAs: 'calls/*')

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('gawk'), eval("awk --version | sed '1!d; s/^GNU Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    methylation_summary.sh -f bedgraph ${args} -- "${prefix}" calls/*
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.methylation_generalstats_mqc.tsv
    touch ${prefix}.methylation_beta_hist_mqc.tsv
    touch ${prefix}.methylation_contigs_mqc.tsv
    """
}
