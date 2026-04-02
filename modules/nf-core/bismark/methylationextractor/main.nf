nextflow.preview.types = true

process BISMARK_METHYLATIONEXTRACTOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38e61d14ccaed82f60c967132963eb467d0fa4bccb7a21404c49b4f377735f03/data' :
        'community.wave.seqera.io/library/bismark:0.25.1--1f50935de5d79c47' }"

    input:
    record(
        meta: Record,
        bam: Path,
        bismark_index: Path
    )

    output:
    record(
        id                   : meta.id,
        meta                 : meta,
        methylation_bedgraph : file("*.bedGraph.gz"),
        methylation_calls    : files("*.txt.gz"),
        methylation_coverage : file("*.cov.gz"),
        methylation_report   : file("*_splitting_report.txt"),
        methylation_mbias    : file("*.M-bias.txt"),
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    // Assign sensible numbers for multicore and buffer_size based on bismark docs
    if(!args.contains('--multicore') && task.cpus >= 6){
        args += " --multicore ${(task.cpus / 3) as int}"
    }
    // Only set buffer_size when there are more than 6.GB of memory available
    if(!args.contains('--buffer_size') && task.memory?.toGiga() > 6){
        args += " --buffer_size ${task.memory.toGiga() - 2}G"
    }

    def seqtype  = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor \\
        ${bam} \\
        --bedGraph \\
        --counts \\
        --gzip \\
        --report \\
        ${seqtype} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.bedGraph.gz
    echo | gzip > ${prefix}.txt.gz
    echo | gzip > ${prefix}.cov.gz
    touch ${prefix}_splitting_report.txt
    touch ${prefix}.M-bias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
