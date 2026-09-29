nextflow.enable.types = true

process PICARD_CREATESEQUENCEDICTIONARY {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data' :
        'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6' }"

    input:
    fasta: Path
    args: String?

    output:
    file("*.dict")

    topic:
    tuple(task.process, 'picard', eval("picard CreateSequenceDictionary --version 2>&1 | sed -n 's/.*Version://p'")) >> 'versions'

    script:
    args = args ?: ''
    def prefix = "${fasta.baseName}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CreateSequenceDictionary] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.toMega()*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CreateSequenceDictionary  \\
        $args \\
        --REFERENCE $fasta \\
        --OUTPUT ${prefix}.dict
    """

    stub:
    def prefix = "${fasta.baseName}"
    """
    touch ${prefix}.dict
    """

}
