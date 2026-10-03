nextflow.enable.types = true

process BWAMETH_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.9--pyh7e72e81_0' :
        'quay.io/biocontainers/bwameth:0.2.9--pyh7e72e81_0' }"

    input:
    record(
        meta: Record,
        reads: List<Path>,
        fasta: Path,
        bwameth_index: Path,
        args: String?,
        args2: String?,
        prefix: String?
    )

    output:
    record(
        meta : meta,
        bam  : file("*.bam")
    )

    topic:
    tuple(task.process, 'bwameth', eval("bwameth.py --version | cut -f2 -d ' '")) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    args = task.ext.args ?: args ?: ''
    args2 = task.ext.args2 ?: args2 ?: ''
    prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    """
    export BWA_METH_SKIP_TIME_CHECKS=1
    ln -sf \$(readlink ${fasta}) ${bwameth_index}/${fasta}

    bwameth.py \\
        ${args} \\
        -t ${task.cpus} \\
        --reference ${bwameth_index}/${fasta} \\
        ${reads.join(' ')} \\
        | samtools view ${args2} -@ ${task.cpus} -bhS -o ${prefix}.bam -
    """

    stub:
    prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
