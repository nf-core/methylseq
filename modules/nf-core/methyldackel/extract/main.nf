nextflow.enable.types = true

process METHYLDACKEL_EXTRACT {
    tag id
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7' :
        'quay.io/biocontainers/methyldackel:0.6.1--he4a0461_7' }"

    input:
    record(
        id: String,
        bam: Path,
        bai: Path,
        fasta: Path,
        fai: Path,
        args: String?
    )

    output:
    record(
        id                     : id,
        methyldackel_bedgraph  : files("*.bedGraph", optional: true),
        methyldackel_methylkit : files("*.methylKit", optional: true)
    )

    topic:
    tuple(task.process, 'methyldackel', eval("MethylDackel --version 2>&1 | cut -f1 -d' '")) >> 'versions'

    script:
    args = args ?: ''
    """
    MethylDackel extract \\
        ${args} \\
        ${fasta} \\
        ${bam}
    """

    stub:
    args = args ?: ''
    def out_extension = args.contains('--methylKit') ? 'methylKit' : 'bedGraph'
    """
    touch ${bam.baseName}_CpG.${out_extension}
    """
}
