nextflow.enable.types = true

process RASTAIR_CALL {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15120636da858ba73a2493281bfa418005f08c0ed09369a837c05f3f9e14a4a6/data' :
        'community.wave.seqera.io/library/rastair:0.8.2--bf70eeab4121509c' }"

    input:
    record(
        meta: Record,
        bam: Path,
        bai: Path,
        fasta: Path,
        fai: Path,
        trim_OT: String,
        trim_OB: String
    )

    output:
    record(
        meta         : meta,
        rastair_call : file("*.rastair_call.txt")
    )

    topic:
    tuple(task.process, 'rastair', eval("rastair --version | sed 's/rastair //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nt_OT_to_trim = meta.trim_OT ?: trim_OT
    def nt_OB_to_trim = meta.trim_OB ?: trim_OB

    """
    rastair call \\
        --threads ${task.cpus} \\
        --nOT ${nt_OT_to_trim} \\
        --nOB ${nt_OB_to_trim} \\
        --fasta-file ${fasta} \\
        ${bam} > ${prefix}.rastair_call.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_call.txt
    """
}
