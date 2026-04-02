nextflow.preview.types = true

process CAT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    record(
        meta: Record,
        reads: List<Path>
    )

    stage:
    stageAs reads, "input*/*"

    output:
    record(
        id: meta.id,
        meta: meta,
        reads: files("*.merged.fastq.gz").toSorted()
    )

    topic:
    file("versions.yml") >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads.collect { file -> "${file}" }.toList()
    if (meta.single_end) {
        if (readList.size() >= 1) {
            """
            cat ${readList.join(' ')} > ${prefix}.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (readList.size() >= 2) {
            def read1 = []
            def read2 = []
            readList.withIndex().each { v, ix -> (ix & 1 ? read2 : read1) << v }
            """
            cat ${read1.join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.join(' ')} > ${prefix}_2.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads.collect { file -> "${file}" }
    if (meta.single_end) {
        if (readList.size() >= 1) {
            """
            echo '' | gzip > ${prefix}.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (readList.size() >= 2) {
            """
            echo '' | gzip > ${prefix}_1.merged.fastq.gz
            echo '' | gzip > ${prefix}_2.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }
}
