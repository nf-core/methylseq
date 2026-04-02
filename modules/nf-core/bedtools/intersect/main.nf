nextflow.preview.types = true

process BEDTOOLS_INTERSECT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    record(
        meta: Record,
        intervals1: Path,
        intervals2: Path,
        chrom_sizes: Path?
    )

    output:
    record(
        id: meta.id,
        meta: meta,
        bedgraph_intersect: file("*.${extension}")
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //Extension of the output file. It is set by the user via "ext.suffix" in the config. Corresponds to the file format which depends on arguments (e. g., ".bed", ".bam", ".txt", etc.).
    extension = task.ext.suffix ?: "${intervals1.extension}"
    def sizes = chrom_sizes ? "-g ${chrom_sizes}" : ''
    if ("$intervals1" == "${prefix}.${extension}" ||
        "$intervals2" == "${prefix}.${extension}")
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    bedtools \\
        intersect \\
        -a $intervals1 \\
        -b $intervals2 \\
        $args \\
        $sizes \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.suffix ?: "bed"
    if ("$intervals1" == "${prefix}.${extension}" ||
        "$intervals2" == "${prefix}.${extension}")
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
