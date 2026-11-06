/*
 * Summarise the rastair version 2 CpG BED as MultiQC custom content:
 * General Statistics columns, a per-CpG methylation distribution and spike-in/control methylation.
 */
process RASTAIR2_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50c87ecbec7ea76f99ff4ad84e16f083c0b74651bf813d9158b1004234f6feca/data'
        : 'community.wave.seqera.io/library/htslib_rastair:4ffd26b40813e268'}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rastair2_summary.sh ${prefix} ${bed} ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair2_generalstats_mqc.tsv
    touch ${prefix}.rastair2_beta_hist_mqc.tsv
    touch ${prefix}.rastair2_contigs_mqc.tsv
    """
}
