nextflow.enable.types = true

process BEDTOOLS_INTERSECT {
    tag id
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    record(
        id: String,
        intervals1: Path,
        intervals2: Path,
        args: String?,
        prefix: String?,
        suffix: String?
    )
    chrom_sizes: Path?

    output:
    record(
        id                 : id,
        coverage_intersect : file("*.${extension}")
    )

    topic:
    tuple(task.process, 'bedtools', eval("bedtools --version | sed -e 's/bedtools v//g'")) >> 'versions'

    script:
    args = args ?: ''
    prefix = prefix ?: "${id}"
    //Extension of the output file. It is set by the user via "ext.suffix" in the config. Corresponds to the file format which depends on arguments (e. g., ".bed", ".bam", ".txt", etc.).
    extension = suffix ?: "${intervals1.extension}"
    def sizes = chrom_sizes ? "-g ${chrom_sizes}" : ''
    if ("$intervals1" == "${prefix}.${extension}" ||
        "$intervals2" == "${prefix}.${extension}")
        error "Input and output names are the same, use \"prefix\" to disambiguate!"
    """
    bedtools \\
        intersect \\
        -a $intervals1 \\
        -b $intervals2 \\
        $args \\
        $sizes \\
        > ${prefix}.${extension}
    """

    stub:
    prefix = prefix ?: "${id}"
    extension = suffix ?: "bed"
    if ("$intervals1" == "${prefix}.${extension}" ||
        "$intervals2" == "${prefix}.${extension}")
        error "Input and output names are the same, use \"prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}
    """
}
