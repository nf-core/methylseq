process TRIMDIVERSITY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=2.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:2.7' :
        'biocontainers/python:2.7' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*_trimmed.log")     , emit: log     , optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        def args_list = args.split("\\s(?=--)").toList()
        args_list.removeAll { it.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz

        trimRRBSdiversityAdaptCustomers.py -1 ${prefix}.fastq.gz

        cp .command.log ${prefix}_trimmed.log

        if [ -f ${prefix}.fastq_trimmed.fq.gz ]; then
            mv ${prefix}.fastq_trimmed.fq.gz ${prefix}_trimmed.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

        trimRRBSdiversityAdaptCustomers.py -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz

        cp .command.log ${prefix}_trimmed.log

        if [ -f ${prefix}_1.fastq_trimmed.fq.gz ]; then
            mv ${prefix}_1.fastq_trimmed.fq.gz ${prefix}_1_trimmed.fastq.gz
        fi
        if [ -f ${prefix}_2.fastq_trimmed.fq.gz ]; then
            mv ${prefix}_2.fastq_trimmed.fq.gz ${prefix}_2_trimmed.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
}
