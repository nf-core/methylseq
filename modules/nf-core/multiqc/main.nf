nextflow.enable.types = true

process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1bef8af6be88c5733461959c46ac8ef73d18f65277f62a1695d0e1633054f9c2/data'
        : 'community.wave.seqera.io/library/multiqc:1.34--db7c73dae76bc9e6'}"

    input:
    record(
        multiqc_files: Set<Path>,
        multiqc_config: List<Path>,
        multiqc_logo: Path?,
        replace_names: Path?,
        sample_names: Path?
    )

    stage:
    stageAs multiqc_files, "?/*"
    stageAs multiqc_config, "?/*"

    output:
    record(
        report : file("*.html"),
        data   : file("*_data"),
        plots  : file("*_plots", optional: true)
    )


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config ${multiqc_config.join(' --config ')}" : ""
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}
