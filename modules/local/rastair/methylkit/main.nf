/**
This module describes the custom Rastair process for
parsing the rastair call output
and converting into MethylKit and Bismark digestible formats.
*/

process CONVERT_TO_METHYLKIT {
    label 'process_low'
    container "docker.io/sbludwig/rastair:version-0.8.2"

    input:
    tuple val(meta), path(rastair_call_txt)

    output:
    tuple val(meta), path("*methylkit.txt.gz"), emit: methylkit
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat ${rastair_call_txt} | /app/scripts/rastair_call_to_methylkit.sh | gzip -c > ${prefix}.rastair_methylkit.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
