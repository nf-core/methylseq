nextflow.enable.types = true

process RASTAIR_MBIASPARSER {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15120636da858ba73a2493281bfa418005f08c0ed09369a837c05f3f9e14a4a6/data' :
        'community.wave.seqera.io/library/rastair:0.8.2--bf70eeab4121509c' }"

    input:
    record(
        meta: Record,
        rastair_mbias: Path
    )

    output:
    record(
        meta              : meta,
        rastair_mbias_pdf : file("*.rastair_mbias_processed.pdf", optional: true),
        rastair_mbias_csv : file("*.rastair_mbias_processed.csv"),
        trim_OT           : env('trim_OT'),
        trim_OB           : env('trim_OB')
    )

    topic:
    tuple(task.process, 'rastair', eval("rastair --version | sed 's/rastair //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_mbias.R --pdf -o ${prefix}.rastair_mbias_processed.pdf ${rastair_mbias} > ${prefix}.rastair_mbias_processed.txt

    parse_mbias.R ${prefix}.rastair_mbias_processed.txt ${prefix}.rastair_mbias_processed.csv
    trim_OT=\$(head -n 1 ${prefix}.rastair_mbias_processed.csv)
    trim_OB=\$(head -n 2 ${prefix}.rastair_mbias_processed.csv | tail -n 1)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_mbias_processed.pdf
    touch ${prefix}.rastair_mbias_processed.csv
    """
}
