// Filter a bedGraph (against widened targets) or a rastair CpG BED (strand rows padded to their CpG) to the targets
process FILTER_BEDGRAPH_TARGETS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(bedgraph), path(targets, stageAs: 'targets/*')

    output:
    tuple val(meta), path("*.targeted.{bedGraph,bed}"), emit: intersect
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (bedgraph.name =~ /\.(bed\.gz|rastair_call\.txt)$/) {
        """
        # Pad each strand row to its CpG (minus-strand rows to the left); keep # header lines, drop track lines
        case "${bedgraph}" in
            *.gz) zcat ${bedgraph} ;;
            *) cat ${bedgraph} ;;
        esac \\
            | awk 'BEGIN{FS=OFS="\t"} /^track/{next} /^#/{print; next} {if(\$6=="-") \$2=\$2-1; else if(NF>=3) \$3=\$3+1; print}' > extended.bedGraph

        # Intersect with target regions (-wa: write original -a entry; -u: unique hits only; -header: keep the -a header)
        bedtools intersect \\
            -header \\
            -a extended.bedGraph \\
            -b ${targets} \\
            -wa \\
            -u \\
            > intersected.bedGraph

        # Restore original single-base coordinates
        awk 'BEGIN{FS=OFS="\t"} /^#/{print; next} {if(\$6=="-") \$2=\$2+1; else if(NF>=3) \$3=\$3-1; print}' intersected.bedGraph \\
            > ${prefix}.targeted.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
    else {
        """
        # Keep # header lines, drop track lines; -wa -u writes each row once, unchanged
        case "${bedgraph}" in
            *.gz) zcat ${bedgraph} ;;
            *) cat ${bedgraph} ;;
        esac \\
            | awk '!/^track/' \\
            | bedtools intersect -header -a stdin -b ${targets} -wa -u \\
            > ${prefix}.targeted.bedGraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
