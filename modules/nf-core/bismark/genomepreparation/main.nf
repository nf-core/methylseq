nextflow.enable.types = true

process BISMARK_GENOMEPREPARATION {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bddea334e6ccbce005ce540214747acf822b040185d2198220dcfbb4b258c331/data' :
        'community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4' }"

    input:
    fasta: Path
    args: String?

    stage:
    stageAs fasta, "BismarkIndex/"

    output:
    file("BismarkIndex")

    topic:
    tuple(task.process, 'bismark', eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'")) >> 'versions'

    script:
    args = args ?: ''
    """
    bismark_genome_preparation \\
        ${args} \\
        BismarkIndex
    """

    stub:
    """
    rm $fasta

    mkdir -p BismarkIndex/Bisulfite_Genome/CT_conversion
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.2.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.3.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.4.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.rev.1.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.rev.2.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa

    mkdir -p BismarkIndex/Bisulfite_Genome/GA_conversion
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.2.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.3.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.4.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.rev.1.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.rev.2.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa
    """
}
