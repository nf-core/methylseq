nextflow.enable.types = true

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    record(
        meta: Record,
        reads: List<Path>,
        bwa_index: Path,
        fasta: Path?
    )
    sort_bam: Boolean

    output:
    record(
        meta : meta,
        bam  : file("*.bam", optional: true),
        cram : file("*.cram", optional: true),
        sam  : file("*.sam", optional: true),
        csi  : file("*.csi", optional: true),
        crai : file("*.crai", optional: true)
    )

    topic:
    tuple(task.process, 'bwa', eval("bwa 2>&1 | sed -n \"s/^Version: //p\"")) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    //
    // For SAM output we can skip samtools view
    //
    def pipe_command = ""
    if (extension == "sam") {
        pipe_command = "> ${prefix}.${extension}"
    } else {
        pipe_command = "| samtools $samtools_command $args2 ${reference} --threads $task.cpus -o ${prefix}.${extension} -"
    }
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        ${reads.join(' ')} \\
        $pipe_command
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.csi
    touch ${prefix}.crai
    """
}
