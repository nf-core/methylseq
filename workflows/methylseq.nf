/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'

// Aligner: bismark or bismark_hisat
if( params.aligner =~ /bismark/ ){
    include { BISMARK } from '../subworkflows/local/bismark'
}
// Aligner: bwameth
else if ( params.aligner == 'bwameth' ){
    include { BWAMETH } from '../subworkflows/local/bwameth'
}

include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC              } from '../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methylseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    /*
     * MODULE: Run TrimGalore!
     */
    if (!params.skip_trimming) {
        TRIMGALORE(ch_samplesheet)
        reads = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    } else {
        reads = ch_samplesheet
    }

    /*
     * SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
     */

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){

        /*
         * Run Bismark alignment + downstream processing
         */
        BISMARK (
            reads,
            PREPARE_GENOME.out.bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq
        )
        ch_versions = ch_versions.mix(BISMARK.out.versions.unique{ it.baseName })
        ch_bam = BISMARK.out.bam
        ch_dedup = BISMARK.out.dedup
        ch_aligner_mqc = BISMARK.out.mqc
    }
    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){

        BWAMETH (
            reads,
            PREPARE_GENOME.out.bwameth_index,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fasta_index,
            params.skip_deduplication || params.rrbs,
        )
        ch_versions = ch_versions.mix(BWAMETH.out.versions.unique{ it.baseName })
        ch_bam = BWAMETH.out.bam
        ch_dedup = BWAMETH.out.dedup
        ch_aligner_mqc = BWAMETH.out.mqc
    }

    /*
     * MODULE: Qualimap BamQC
     */
    QUALIMAP_BAMQC (
        ch_dedup,
        params.bamqc_regions_file ? Channel.fromPath( params.bamqc_regions_file, checkIfExists: true ).toList() : []
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    /*
     * MODULE: Run Preseq
     */
    PRESEQ_LCEXTRAP (
        ch_bam
    )
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{ it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect{ it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc.ifEmpty([]))
    if (!params.skip_trimming) {
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ it[1] })
    }
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
