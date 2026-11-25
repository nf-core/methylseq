// Merge several target BEDs into their union, for tools that take a single target set (riker hybcap)
process MERGE_TARGETS {
    tag "${beds.size()} panels"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    path beds, stageAs: 'targets/?/*'

    output:
    path "merged_targets.bed", emit: targets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # The image's gzip has no -f for plain files; drop track, browser and # lines and any extra columns
    for bed in ${beds.join(' ')}; do if gzip -cd "\$bed" > /dev/null 2>&1; then gzip -cd "\$bed"; else cat "\$bed"; fi; done \\
        | tr -d '\\r' \\
        | awk 'BEGIN{FS=OFS="\\t"} /^(track|browser|#)/ || NF < 3 {next} {print \$1, \$2, \$3}' \\
        | LC_ALL=C sort -k1,1 -k2,2n \\
        | bedtools merge -i stdin \\
        > merged_targets.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    """
    touch merged_targets.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
