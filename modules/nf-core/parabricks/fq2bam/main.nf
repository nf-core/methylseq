nextflow.enable.types = true

process PARABRICKS_FQ2BAM {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to run on a cluster because we need to copy the fasta reference, see https://github.com/nf-core/modules/issues/9230
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    record(
        meta: Record,
        reads: List<Path>,
        fasta: Path,
        bwa_index: Path
    )
    intervals: List<Path>
    known_sites: List<Path>
    output_fmt: String

    output:
    record(
        meta              : meta,
        bam               : file("*.bam", optional: true),
        bai               : file("*.bai", optional: true),
        cram              : file("*.cram", optional: true),
        crai              : file("*.crai", optional: true),
        bqsr_table        : file("*.table", optional: true),
        qc_metrics        : file("*_qc_metrics", optional: true),
        duplicate_metrics : file("*.duplicate-metrics.txt", optional: true)
    )

    topic:
    tuple(task.process, 'parabricks', eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'")) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def in_fq_command = meta.single_end ? "--in-se-fq ${reads.join(' ')}" : "--in-fq ${reads.join(' ')}"
    def extension = "${output_fmt}"

    def known_sites_command = known_sites.collect { knownSite -> "--knownSites ${knownSite}" }.join(' ')
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def intervals_command = intervals.collect { interval -> "--interval-file ${interval}" }.join(' ')

    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp ${fasta} \$INDEX

    pbrun \\
        fq2bam \\
        --ref \$INDEX \\
        ${in_fq_command} \\
        --out-bam ${prefix}.${extension} \\
        ${known_sites_command} \\
        ${known_sites_output_cmd} \\
        ${intervals_command} \\
        ${num_gpus} \\
        --bwa-cpu-thread-pool ${task.cpus} \\
        --monitor-usage \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "${output_fmt}"
    def extension_index = "${output_fmt}" == "cram" ? "crai" : "bai"
    def known_sites_output = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.${extension}
    touch ${prefix}.${extension}.${extension_index}
    ${known_sites_output}
    ${qc_metrics_output}
    ${duplicate_metrics_output}

    # Capture once and build single-line compatible_with (spaces only, no tabs)
    pbrun_version_output=\$(pbrun fq2bam --version 2>&1)

    # Because of a space between BWA and mem in the version output this is handled different to the other modules
    compat_line=\$(echo "\$pbrun_version_output" | awk -F':' '
        /Compatible With:/ {on=1; next}
        /^---/ {on=0}
        on && /:/ {
            key=\$1; val=\$2
            gsub(/[ \\t]+/, " ", key); gsub(/^[ \\t]+|[ \\t]+\$/, "", key)
            gsub(/[ \\t]+/, " ", val); gsub(/^[ \\t]+|[ \\t]+\$/, "", val)
            a[++i]=key ": " val
        }
        END { for (j=1;j<=i;j++) printf "%s%s", (j>1?", ":""), a[j] }
    ')

    cat <<EOF > compatible_versions.yml
    "${task.process}":
    pbrun_version: \$(echo "\$pbrun_version_output" | awk '/^pbrun:/ {print \$2; exit}')
    compatible_with: "\$compat_line"
    EOF
    """
}
