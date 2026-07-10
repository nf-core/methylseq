/*
 * Filter bedGraph files against target regions with CpG-aware boundary handling.
 *
 * Bismark bedGraph files use 0-based half-open coordinates for single-C positions
 * (e.g., chr1 10788 10789 for a CpG at 1-based position 10789). When a CpG
 * straddles a target boundary — C just outside, G just inside — a naive
 * intersection drops it because the single-C interval [10788, 10789) does not
 * overlap the target [10789, …).
 *
 * This process extends each bedGraph entry by 1 bp on the right to represent the
 * full CpG dinucleotide, intersects with the target BED, then restores the
 * original single-base coordinates.
 */
process FILTER_BEDGRAPH_TARGETS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(bedgraph), path(targets)

    output:
    tuple val(meta), path("*.targeted.bedGraph"), emit: intersect
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Read the bedGraph (gzipped from Bismark, plain text from MethylDackel) and extend the end
    # coordinate by 1 bp so the interval covers the full CpG dinucleotide.
    # Skip any track/header lines to avoid corrupting them
    case "${bedgraph}" in
        *.gz) zcat ${bedgraph} ;;
        *) cat ${bedgraph} ;;
    esac \\
        | awk 'BEGIN{OFS="\t"} /^track/{print; next} {if(NF>=3) \$3=\$3+1; print}' > extended.bedGraph

    # Intersect with target regions (-wa: write original -a entry; -u: unique hits only)
    bedtools intersect \\
        -a extended.bedGraph \\
        -b ${targets} \\
        -wa \\
        -u \\
        > intersected.bedGraph

    # Restore original single-base end coordinates
    awk 'BEGIN{OFS="\t"} /^track/{print; next} {if(NF>=3) \$3=\$3-1; print}' intersected.bedGraph \\
        > ${prefix}.targeted.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
