nextflow.enable.types = true

process TRIMGALORE {
    tag id
    label 'process_medium'
    label 'process_low_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e00369598bd6b7b34a7c83d5496c381104bf8b885c31a4b65b92e6ea2059fbb3/data' :
        'community.wave.seqera.io/library/trim-galore:2.3.0--6a38a479b4972363'}"

    input:
    record(
        id: String,
        single_end: Boolean,
        reads: List<Path>,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id            : id,
        trim_reads    : files("*{3prime,5prime,trimmed,val}{,_1,_2}.fq.gz").toSorted(),
        trim_log      : files("*report.txt", optional: true).toSorted(),
        trim_unpaired : files("*unpaired{,_1,_2}.fq.gz", optional: true).toSorted(),
        trim_html     : files("*.html", optional: true).toSorted(),
        trim_zip      : files("*.zip", optional: true).toSorted()
    )

    topic:
    tuple(task.process, 'trimgalore', eval('trim_galore --version | grep -Eo "[0-9]+(\\.[0-9]+)+"')) >> 'versions'

    script:
    args = args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (single_end) {
            cores = (task.cpus as int) - 3
        }
        if (cores < 1) {
            cores = 1
        }
        if (cores > 8) {
            cores = 8
        }
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    prefix = prefix ?: "${id}"
    if (single_end) {
        def args_se = args.replaceAll(/(?i)--\S*_r2\s+\S+/, '').trim()
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s ${reads[0]} ${prefix}.fastq.gz
        trim_galore \\
            ${args_se} \\
            --cores ${cores} \\
            --gzip \\
            ${prefix}.fastq.gz
        """
    }
    else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trim_galore \\
            ${args} \\
            --cores ${cores} \\
            --paired \\
            --gzip \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz
        """
    }

    stub:
    prefix = prefix ?: "${id}"
    if (single_end) {
        output_command = "echo '' | gzip > ${prefix}_trimmed.fq.gz ;"
        output_command += "touch ${prefix}.fastq.gz_trimming_report.txt"
    }
    else {
        output_command = "echo '' | gzip > ${prefix}_1_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_1.fastq.gz_trimming_report.txt ;"
        output_command += "echo '' | gzip > ${prefix}_2_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_2.fastq.gz_trimming_report.txt"
    }
    """
    ${output_command}
    """
}
