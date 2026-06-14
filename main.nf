#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/methylseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/methylseq
    Website: https://nf-co.re/methylseq
    Slack  : https://nfcore.slack.com/channels/methylseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_INDEX_METHYLSEQ     } from './subworkflows/nf-core/fasta_index_methylseq/main'
include { BWA_INDEX                 } from './modules/nf-core/bwa/index/main'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { getGenomeAttribute        } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { METHYLSEQ                 } from './workflows/methylseq/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta         = getGenomeAttribute('fasta')
params.fasta_index   = getGenomeAttribute('fasta_index')
params.bwameth_index = getGenomeAttribute('bwameth')
params.bwamem_index  = getGenomeAttribute('bwa')
params.bismark_index = params.aligner == 'bismark_hisat' ? getGenomeAttribute('bismark_hisat2') : getGenomeAttribute('bismark')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_METHYLSEQ {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = channel.empty()

    //
    // Initialize file channels or values based on params
    //
    ch_fasta                = params.fasta         ? channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] } : channel.empty()
    ch_or_val_fasta_index   = params.fasta_index   ? channel.fromPath(params.fasta_index).map{ it -> [ [id:it.baseName], it ] } : []
    ch_or_val_bismark_index = params.bismark_index ? channel.fromPath(params.bismark_index).map{ it -> [ [id:it.baseName], it ] } : []
    ch_or_val_bwameth_index = params.bwameth_index ? channel.fromPath(params.bwameth_index).map{ it -> [ [id:it.baseName], it ] } : []
    ch_or_val_bwamem_index  = params.bwamem_index  ? channel.fromPath(params.bwamem_index).map{ it -> [ [id:it.baseName], it ] } : []

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    FASTA_INDEX_METHYLSEQ(
        ch_fasta,
        ch_or_val_fasta_index,
        ch_or_val_bismark_index,
        ch_or_val_bwameth_index,
        ch_or_val_bwamem_index,
        params.aligner,
        params.collecthsmetrics,
        params.use_mem2
    )
    ch_versions = ch_versions.mix(FASTA_INDEX_METHYLSEQ.out.versions)

    //
    // WORKFLOW: Run pipeline
    //

    METHYLSEQ (
        samplesheet,
        ch_versions,
        FASTA_INDEX_METHYLSEQ.out.fasta,
        FASTA_INDEX_METHYLSEQ.out.fasta_index,
        FASTA_INDEX_METHYLSEQ.out.bismark_index,
        FASTA_INDEX_METHYLSEQ.out.bwameth_index,
        FASTA_INDEX_METHYLSEQ.out.bwamem_index,
    )
    ch_versions = ch_versions.mix(METHYLSEQ.out.versions)

    emit:
    multiqc_report = METHYLSEQ.out.multiqc_report // channel: [ path(multiqc_report.html )  ]
    versions       = ch_versions                  // channel: [ path(versions.yml) ]

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_METHYLSEQ (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_METHYLSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
