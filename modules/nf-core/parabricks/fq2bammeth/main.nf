nextflow.preview.types = true

process PARABRICKS_FQ2BAMMETH {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.2-1"

    input:
    record(
        meta: Record,
        reads: Path,
        fasta: Path,
        bwameth_index: Path,
        known_sites: Path?
    )

    output:
    record(
        id                : meta.id,
        meta              : meta,
        bam               : file("*.bam"),
        bai               : file("*.bai"),
        qc_metrics        : file("qc_metrics", optional: true),
        bqsr_table        : file("*.table", optional: true),
        duplicate_metrics : file("duplicate-metrics.txt", optional: true),
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args                = task.ext.args ?: ''
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def in_fq_command       = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def known_sites_command = known_sites ? "--knownSites ${known_sites}" : ""
    def known_sites_output  = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def num_gpus            = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    if [ -L $fasta ]; then
        ln -sf \$(readlink $fasta) ${bwameth_index}/$fasta
    else
        ln -sf ../$fasta ${bwameth_index}/$fasta
    fi

    pbrun \\
        fq2bam_meth \\
        --ref ${bwameth_index}/$fasta \\
        $in_fq_command \\
        --out-bam ${prefix}.bam \\
        $known_sites_command \\
        $known_sites_output \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
