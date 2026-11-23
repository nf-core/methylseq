// Widen each target by 1 bp at an edge that cuts a reference CpG, so filtered calls keep both strands of it
process WIDEN_TARGETS_CPG {
    tag "${targets.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(targets, stageAs: 'targets/*')
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("widened_targets.bed"), emit: targets
    tuple val(meta), path("target_issues.txt"), emit: issues
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    widen_targets_cpg.sh ${targets} ${fasta} ${fai}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    """
    touch widened_targets.bed target_issues.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
