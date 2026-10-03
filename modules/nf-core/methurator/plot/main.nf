nextflow.enable.types = true

process METHURATOR_PLOT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/methurator:2.2.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/methurator:2.2.0--pyhdfd78af_0'}"

    input:
    record(
        meta: Record,
        methurator_summary: Path
    )

    output:
    record(
        meta             : meta,
        methurator_plots : files("plots/*.html")
    )

    topic:
    tuple(task.process, 'methurator', eval("methurator --version | sed 's/.* //'")) >> 'versions'

    script:
    """
    methurator plot \\
        --summary ${methurator_summary} \\
        --outdir .

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir plots/
    touch plots/${prefix}.html

    """
}
