/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { FASTQC                    } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE                } from '../../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC            } from '../../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP           } from '../../modules/nf-core/preseq/lcextrap/main'
include { MULTIQC                   } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                 } from '../../modules/nf-core/cat/fastq/main'
include { FASTQ_ALIGN_DEDUP_BISMARK } from '../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { FASTQ_ALIGN_DEDUP_BWAMETH } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth/main'
include { FASTQ_ALIGN_DEDUP_BWAMEM  } from '../../subworkflows/local/fastq_align_dedup_bwamem/main'
include { paramsSummaryMultiqc      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { validateInputSamplesheet  } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { TAPS_CONVERSION           } from '../../subworkflows/local/taps_conversion'
include { METHYLDACKEL              } from '../../subworkflows/local/methyldackel/main'
include { TARGETED_SEQUENCING       } from '../../subworkflows/local/targeted_sequencing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {

    take:
    samplesheet        // channel: [ path(samplesheet.csv) ]
    ch_versions        // channel: [ path(versions.yml)    ]
    ch_fasta           // channel: [ path(fasta)           ]
    ch_fasta_index     // channel: [ path(fasta index)     ]
    ch_bismark_index   // channel: [ path(bismark index)   ]
    ch_bwameth_index   // channel: [ path(bwameth index)   ]
    ch_bwamem_index    // channel: [ path(bwamem_index)    ]

    main:
    ch_fastq         = Channel.empty()
    ch_fastqc_html   = Channel.empty()
    ch_fastqc_zip    = Channel.empty()
    ch_reads         = Channel.empty()
    ch_bam           = Channel.empty()
    ch_bai           = Channel.empty()
    ch_bedgraph      = Channel.empty()
    ch_aligner_mqc   = Channel.empty()
    ch_qualimap      = Channel.empty()
    ch_preseq        = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Branch channels from input samplesheet channel
    //
    ch_samplesheet = samplesheet
                        .branch { meta, fastqs ->
                            single  : fastqs.size() == 1
                                return [ meta, fastqs.flatten() ]
                            multiple: fastqs.size() > 1
                                return [ meta, fastqs.flatten() ]
                        }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_samplesheet.multiple
    )
    ch_fastq    = CAT_FASTQ.out.reads.mix(ch_samplesheet.single)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC (
            ch_fastq
        )
        ch_fastqc_html   = FASTQC.out.html
        ch_fastqc_zip    = FASTQC.out.zip
        ch_versions      = ch_versions.mix(FASTQC.out.versions)
    } else {
        ch_fastqc_html   = Channel.empty()
        ch_fastqc_zip    = Channel.empty()
    }

    //
    // MODULE: Run TrimGalore!
    //
    if (!params.skip_trimming) {
        TRIMGALORE(
            ch_fastq
        )
        ch_reads    = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    } else {
        ch_reads    = ch_fastq
    }

    //
    // SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
    //

    if (params.taps & params.aligner != 'bwamem') {
        log.info "TAPS protocol detected and aligner is not 'bwamem'. We recommend using bwa-mem for TAPS protocol as it is optimized for this type of data."
        // params.aligner = 'bwamem'
    } 

    // Aligner: bismark or bismark_hisat
    if (params.aligner =~ /bismark/ ) {
        //
        // Run Bismark alignment + downstream processing
        //
        ch_bismark_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_bismark_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_bismark, bismark_index ->
                reads: [ meta, reads ]
                fasta: [ meta_fasta, fasta ]
                bismark_index: [ meta_bismark, bismark_index ]
            }

        FASTQ_ALIGN_DEDUP_BISMARK (
            ch_bismark_inputs.reads,
            ch_bismark_inputs.fasta,
            ch_bismark_inputs.bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq
        )
        ch_bam         = FASTQ_ALIGN_DEDUP_BISMARK.out.bam
        ch_bai         = FASTQ_ALIGN_DEDUP_BISMARK.out.bai
        ch_bedgraph    = FASTQ_ALIGN_DEDUP_BISMARK.out.methylation_bedgraph
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BISMARK.out.multiqc
        ch_versions    = ch_versions.mix(FASTQ_ALIGN_DEDUP_BISMARK.out.versions)
    }
    // Aligner: bwameth
    else if (params.aligner == 'bwameth'){

        ch_bwameth_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_bwameth_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_fasta_index, fasta_index, meta_bwameth, bwameth_index ->
                reads: [ meta, reads ]
                fasta: [ meta_fasta, fasta ]
                fasta_index: [ meta_fasta_index, fasta_index ]
                bwameth_index: [ meta_bwameth, bwameth_index ]
            }

        FASTQ_ALIGN_DEDUP_BWAMETH (
            ch_bwameth_inputs.reads,
            ch_bwameth_inputs.fasta,
            ch_bwameth_inputs.fasta_index,
            ch_bwameth_inputs.bwameth_index,
            params.skip_deduplication || params.rrbs,
            workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1
        )
        ch_bam         = FASTQ_ALIGN_DEDUP_BWAMETH.out.bam
        ch_bai         = FASTQ_ALIGN_DEDUP_BWAMETH.out.bai
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BWAMETH.out.multiqc
        ch_versions    = ch_versions.mix(FASTQ_ALIGN_DEDUP_BWAMETH.out.versions)
    }

    // Aligner: bwamem
    else if (params.aligner == 'bwamem'){

        ch_bwamem_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_bwamem_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_fasta_index, fasta_index, meta_bwamem, bwamem_index ->
                reads: [ meta, reads ]
                fasta: [ meta_fasta, fasta ]
                fasta_index: [ meta_fasta_index, fasta_index ]
                bwamem_index: [ meta_bwamem, bwamem_index ]
            }

        FASTQ_ALIGN_DEDUP_BWAMEM (
            ch_bwamem_inputs.reads,
            ch_bwamem_inputs.fasta,
            ch_bwamem_inputs.fasta_index,
            ch_bwamem_inputs.bwamem_index,
            params.skip_deduplication,
            workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1
        )

        ch_bam         = FASTQ_ALIGN_DEDUP_BWAMEM.out.bam
        ch_bai         = FASTQ_ALIGN_DEDUP_BWAMEM.out.bai
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BWAMEM.out.multiqc
        ch_versions    = ch_versions.mix(FASTQ_ALIGN_DEDUP_BWAMEM.out.versions.unique{ it.baseName })
    }   

    else {
        error "ERROR: Invalid aligner '${params.aligner}'. Valid options are: 'bismark', 'bismark_hisat', 'bwameth' or 'bwamem'."
    }

    //
    // Subworkflow: Count positive mC->T conversion rates as a readout for DNA methylation
    //
    if (params.taps || params.aligner == 'bwamem') {
        log.info "TAPS protocol detected. Running TAPS conversion module."
        TAPS_CONVERSION (
            ch_bam,
            ch_bai,
            ch_fasta,
            ch_fasta_index,
        )
        ch_rastair_mbias = TAPS_CONVERSION.out.mbias // channel: [ val(meta), [ txt ] ]
        ch_rastair_call  = TAPS_CONVERSION.out.call // channel: [ val(meta), [ txt ] ]
        ch_versions      = ch_versions.mix(TAPS_CONVERSION.out.versions)
    } 

    //
    // Subworkflow: Count negative C->T conversion rates as a readout for DNA methylation
    //
    else if (!params.taps && (params.aligner == 'bwameth' || (params.aligner == 'bwamem'))) {
        METHYLDACKEL (
            ch_bam,
            ch_bai,
            ch_fasta,
            ch_fasta_index
        )
        ch_bedgraph    = METHYLDACKEL.out.methydackel_extract_bedgraph  // channel: [ val(meta), [ bedgraph ]  ]
        ch_methylkit   = METHYLDACKEL.out.methydackel_extract_methylkit // channel: [ val(meta), [ methylkit ] ]
        ch_mbias       = METHYLDACKEL.out.methydackel_mbias // channel: [ val(meta), [ mbias ] ]
        ch_versions    = ch_versions.mix(METHYLDACKEL.out.versions)
    }

    //
    // MODULE: Qualimap BamQC
    // skipped by default. to use run with `--run_qualimap` param.
    //
    if(params.run_qualimap) {
        QUALIMAP_BAMQC (
            ch_bam,
            params.bamqc_regions_file ? Channel.fromPath( params.bamqc_regions_file, checkIfExists: true ).toList() : []
        )
        ch_qualimap = QUALIMAP_BAMQC.out.results
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
    }

    //
    // MODULE: Targeted sequencing analysis
    // skipped by default. to use run with `--run_targeted_sequencing` param.
    //
    if (params.run_targeted_sequencing){
        if (!params.target_regions_file) {
            error "ERROR: --target_regions_file must be specified when using --run_targeted_sequencing"
        }
        TARGETED_SEQUENCING (
            ch_bedgraph,
            Channel.fromPath(params.target_regions_file, checkIfExists: true),
            ch_fasta,
            ch_fasta_index,
            ch_bam,
            ch_bai,
            params.collecthsmetrics
        )
        ch_versions = ch_versions.mix(TARGETED_SEQUENCING.out.versions)
    }

    //
    // MODULE: Preseq LCEXTRAP
    // skipped by default. to use run with `--run_preseq` param.
    //
    if(params.run_preseq) {
        PRESEQ_LCEXTRAP (
            ch_bam
        )
        ch_preseq   = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions)
    }

    //
    // Collate and save software versions
    //
    ch_collated_versions = softwareVersionsToYAML(ch_versions)
                                .collectFile(
                                    storeDir: "${params.outdir}/pipeline_info",
                                    name: 'nf_core_'  +  'methylseq_software_'  + 'mqc_'  + 'versions.yml',
                                    sort: true,
                                    newLine: true
                                )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            Channel.fromPath(params.multiqc_config, checkIfExists: true) :
            Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        if(params.run_qualimap) {
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{ it[1] }.ifEmpty([]))
        }
        if (params.run_preseq) {
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect{ it[1] }.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc.ifEmpty([]))
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ it[1] })
        }
        if (params.run_targeted_sequencing) {
            if (params.collecthsmetrics) {
                ch_multiqc_files = ch_multiqc_files.mix(TARGETED_SEQUENCING.out.picard_metrics.collect{ it[1] }.ifEmpty([]))
            }
        }
        if (!params.skip_fastqc) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))
        }

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        ch_multiqc_report = MULTIQC.out.report.toList()
    } else {
        ch_multiqc_report = Channel.empty()
    }

    emit:
    bam            = ch_bam                      // channel: [ val(meta), path(bam) ]
    bai            = ch_bai                      // channel: [ val(meta), path(bai) ]
    qualimap       = ch_qualimap                 // channel: [ val(meta), path(qualimap) ]
    preseq         = ch_preseq                   // channel: [ val(meta), path(preseq) ]
    multiqc_report = ch_multiqc_report            // channel: [ path(multiqc_report.html )  ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
