nextflow.preview.types = true

process BISMARK_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38e61d14ccaed82f60c967132963eb467d0fa4bccb7a21404c49b4f377735f03/data' :
        'community.wave.seqera.io/library/bismark:0.25.1--1f50935de5d79c47' }"

    input:
    record(
        meta: Record,
        reads: List<Path>,
        fasta: Path,
        bismark_index: Path
    )

    stage:
    stageAs fasta, 'tmp/*' // This change mounts as directory containing the FASTA file to prevent nested symlinks

    output:
    record(
        id: meta.id,
        meta: meta,
        bam: file("*bam"),
        align_report: file("*report.txt"),
        unmapped: file("*fq.gz", optional: true)
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def args = task.ext.args ?: ''
    if(task.ext.prefix){
        args += " --prefix ${task.ext.prefix}"
    }
    def fastq = meta.single_end ? "${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"

    // Try to assign sensible bismark --multicore if not already set
    if(!args.contains('--multicore') && task.cpus){

        // Numbers based on recommendation by Felix for a typical mouse genome
        def ccore = 1
        def cpu_per_multicore = 3
        def mem_per_multicore = (13.GB).toBytes()
        if(args.contains('--non_directional')){
            cpu_per_multicore = 5
            mem_per_multicore = (18.GB).toBytes()
        }

        // How many multicore splits can we afford with the cpus we have?
        ccore = ((task.cpus as int) / cpu_per_multicore) as int

        // Check that we have enough memory
        try {
            def tmem = (task.memory as MemoryUnit).toBytes()
            def mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (all) {
            log.warn "Not able to define bismark align multicore based on available memory"
        }
        if(ccore > 1){
            args += " --multicore ${ccore}"
        }
    }
    """
    bismark \\
        ${fastq} \\
        --genome ${bismark_index} \\
        --bam \\
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
    touch ${prefix}.bam
    touch ${prefix}.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
