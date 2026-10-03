nextflow.enable.types = true

process BWAMETH_INDEX {
    tag "${fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bwameth:0.2.9--pyh7e72e81_0'
        : 'quay.io/biocontainers/bwameth:0.2.9--pyh7e72e81_0'}"

    input:
    fasta: Path
    use_mem2: Boolean

    stage:
    stageAs fasta, "BwamethIndex/"

    output:
    file("BwamethIndex")

    topic:
    tuple(task.process, 'bwameth', eval("bwameth.py --version | cut -f2 -d' '")) >> 'versions'

    script:
    def index_cmd = use_mem2 ? "index-mem2" : "index"
    """
    bwameth.py ${index_cmd} ${fasta}

    rm $fasta
    """

    stub:
    """
    rm ${fasta}

    mkdir -p BwamethIndex/
    touch BwamethIndex/genome.fasta.bwameth.c2t
    touch BwamethIndex/genome.fasta.bwameth.c2t.amb
    touch BwamethIndex/genome.fasta.bwameth.c2t.ann
    touch BwamethIndex/genome.fasta.bwameth.c2t.bwt
    touch BwamethIndex/genome.fasta.bwameth.c2t.pac
    touch BwamethIndex/genome.fasta.bwameth.c2t.sa
    """
}
