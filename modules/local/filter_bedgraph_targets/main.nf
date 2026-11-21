/*
 * Filter bedGraph files against target regions with CpG-aware boundary handling.
 *
 * Bismark bedGraph files use 0-based half-open coordinates for single-C positions
 * (e.g., chr1 10788 10789 for a CpG at 1-based position 10789). When a CpG
 * straddles a target boundary — C just outside, G just inside — a naive
 * intersection drops it because the single-C interval [10788, 10789) does not
 * overlap the target [10789, …).
 *
 * This process extends each entry by 1 bp to represent the full CpG dinucleotide,
 * intersects with the target BED, then restores the original single-base coordinates.
 * Entries extend to the right, except minus-strand rows of a rastair 2 CpG BED
 * (the G of each CpG), which extend to the left. A rastair 2 BED is written as
 * *.targeted.bed with its # header; bedGraphs as *.targeted.bedGraph.
 */
process FILTER_BEDGRAPH_TARGETS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(bedgraph), path(targets)

    output:
    tuple val(meta), path("*.targeted.{bedGraph,bed}"), emit: intersect
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = bedgraph.name.endsWith('.bed.gz') ? 'bed' : 'bedGraph'
    """
    # Read the bedGraph (gzipped from Bismark, plain text from MethylDackel) or rastair 2 CpG BED (gzipped)
    # and extend each single-base entry by 1 bp so it covers the full CpG dinucleotide. rastair 2 writes the
    # G of each CpG as its own minus-strand row (column 6), which extends to the left.
    # Drop track lines and keep # header lines: rastair 2 names its columns in a # header
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
        > ${prefix}.targeted.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
