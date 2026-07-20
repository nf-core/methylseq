/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap              } from 'plugin/nf-schema'
include { FASTQC                        } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE                    } from '../../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC                } from '../../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP               } from '../../modules/nf-core/preseq/lcextrap/main'
include { MULTIQC                       } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                     } from '../../modules/nf-core/cat/fastq/main'
include { FASTQ_ALIGN_DEDUP_BISMARK     } from '../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { FASTQ_ALIGN_DEDUP_BWAMETH     } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth/main'
include { FASTQ_ALIGN_DEDUP_BWAMEM      } from '../../subworkflows/nf-core/fastq_align_dedup_bwamem/main'
include { PICARD_MARKDUPLICATES         } from '../../modules/nf-core/picard/markduplicates/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { paramsSummaryMultiqc          } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { validateInputSamplesheet      } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { BAM_TAPS_CONVERSION           } from '../../subworkflows/nf-core/bam_taps_conversion'
include { BAM_METHYLDACKEL              } from '../../subworkflows/nf-core/bam_methyldackel/main'
include { TARGETED_SEQUENCING           } from '../../subworkflows/local/targeted_sequencing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {
    take:
    samplesheet // channel: [ path(samplesheet.csv) ]
    ch_versions // channel: [ path(versions.yml)    ]
    ch_fasta // channel: [ path(fasta)           ]
    ch_fasta_index // channel: [ path(fasta index)     ]
    ch_bismark_index // channel: [ path(bismark index)   ]
    ch_bwameth_index // channel: [ path(bwameth index)   ]
    ch_bwamem_index // channel: [ path(bwamem_index)    ]

    main:
    ch_fastq = channel.empty()
    ch_fastqc_html = channel.empty()
    ch_fastqc_zip = channel.empty()
    ch_reads = channel.empty()
    ch_bam = channel.empty()
    ch_bai = channel.empty()
    ch_gzi = channel.empty()
    ch_bedgraph = channel.empty()
    ch_coverage = channel.empty()
    ch_aligner_mqc = channel.empty()
    ch_rastair_mbias = channel.empty()
    ch_rastair_call = channel.empty()
    ch_methylkit = channel.empty()
    ch_mbias = channel.empty()
    ch_qualimap = channel.empty()
    ch_preseq = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Branch channels from input samplesheet channel
    //
    ch_samplesheet = samplesheet.branch { meta, fastqs ->
        single: fastqs.size() == 1
        return [meta, fastqs.flatten()]
        multiple: fastqs.size() > 1
        return [meta, fastqs.flatten()]
    }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(
        ch_samplesheet.multiple
    )
    ch_fastq = CAT_FASTQ.out.reads.mix(ch_samplesheet.single)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC(
            ch_fastq
        )
        ch_fastqc_html = FASTQC.out.html
        ch_fastqc_zip = FASTQC.out.zip
    }
    else {
        ch_fastqc_html = channel.empty()
        ch_fastqc_zip = channel.empty()
    }

    //
    // MODULE: Run TrimGalore!
    //
    if (!params.skip_trimming) {
        TRIMGALORE(
            ch_fastq
        )
        ch_reads = TRIMGALORE.out.reads
    }
    else {
        ch_reads = ch_fastq
    }

    //
    // SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
    //

    if (params.taps && params.aligner != 'bwamem') {
        log.info("TAPS protocol detected and aligner is not 'bwamem'. We recommend using bwa-mem for TAPS protocol as it is optimized for this type of data.")
    }

    // Aligner: bismark or bismark_hisat
    if (params.aligner =~ /bismark/) {
        //
        // Run Bismark alignment + downstream processing
        //
        ch_bismark_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_bismark_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_bismark, bismark_index ->
                reads: [meta, reads]
                fasta_fai: [meta_fasta, fasta, []]
                bismark_index: [meta_bismark, bismark_index]
            }

        FASTQ_ALIGN_DEDUP_BISMARK(
            ch_bismark_inputs.reads,
            ch_bismark_inputs.fasta_fai,
            ch_bismark_inputs.bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq,
        )
        ch_bam = FASTQ_ALIGN_DEDUP_BISMARK.out.bam
        ch_bai = FASTQ_ALIGN_DEDUP_BISMARK.out.index
        ch_bedgraph = FASTQ_ALIGN_DEDUP_BISMARK.out.methylation_bedgraph
        ch_coverage = FASTQ_ALIGN_DEDUP_BISMARK.out.methylation_coverage
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BISMARK.out.multiqc
    }
    else if (params.aligner == 'bwameth') {

        ch_bwameth_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_bwameth_index)
            .multiMap { meta, reads, meta_fasta, fasta, _meta_fasta_index, fasta_index, meta_bwameth, bwameth_index ->
                reads: [meta, reads]
                fasta_fai: [meta_fasta, fasta, fasta_index]
                bwameth_index: [meta_bwameth, bwameth_index]
            }

        FASTQ_ALIGN_DEDUP_BWAMETH(
            ch_bwameth_inputs.reads,
            ch_bwameth_inputs.fasta_fai,
            ch_bwameth_inputs.bwameth_index,
            params.skip_deduplication || params.rrbs,
            workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1,
        )
        ch_bam = FASTQ_ALIGN_DEDUP_BWAMETH.out.bam
        ch_bai = FASTQ_ALIGN_DEDUP_BWAMETH.out.bai
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BWAMETH.out.multiqc
    }
    else if (params.aligner == 'bwamem') {

        ch_bwamem_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_bwamem_index)
            .multiMap { meta, reads, meta_fasta, fasta, _meta_fasta_index, fasta_index, meta_bwamem, bwamem_index ->
                reads: [meta, reads]
                fasta_fai: [meta_fasta, fasta, fasta_index]
                bwamem_index: [meta_bwamem, bwamem_index]
            }

        FASTQ_ALIGN_DEDUP_BWAMEM(
            ch_bwamem_inputs.reads,
            ch_bwamem_inputs.fasta_fai,
            ch_bwamem_inputs.bwamem_index,
            params.skip_deduplication,
            workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1,
            'bam',
            [[:], []],
            [[:], []],
        )

        ch_bam = FASTQ_ALIGN_DEDUP_BWAMEM.out.bam
        ch_bai = FASTQ_ALIGN_DEDUP_BWAMEM.out.index
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BWAMEM.out.multiqc
    }
    else {
        error("ERROR: Invalid aligner '${params.aligner}'. Valid options are: 'bismark', 'bismark_hisat', 'bwameth' or 'bwamem'.")
    }

    //
    // Subworkflow: Count positive mC->T conversion rates as a readout for DNA methylation
    //
    if (params.taps || params.aligner == 'bwamem') {

        ch_bam_bai = ch_bam.join(ch_bai)
        ch_taps_inputs = ch_bam_bai
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .multiMap { meta, bam, bai, _meta_fasta, fasta, _meta_fai, fai ->
                bam: [meta, bam]
                bai: [meta, bai]
                fasta: [meta, fasta]
                fasta_index: [meta, fai]
            }

        BAM_TAPS_CONVERSION(
            ch_taps_inputs.bam,
            ch_taps_inputs.bai,
            ch_taps_inputs.fasta,
            ch_taps_inputs.fasta_index,
        )
        ch_rastair_mbias = BAM_TAPS_CONVERSION.out.mbias
        // channel: [ val(meta), [ txt ] ]
        ch_rastair_call = BAM_TAPS_CONVERSION.out.call
        // channel: [ val(meta), [ txt ] ]
        ch_versions = ch_versions.mix(BAM_TAPS_CONVERSION.out.versions)
    }
    else if (!params.taps && params.aligner == 'bwameth') {

        ch_bam_bai = ch_bam.join(ch_bai)
        ch_methyldackel_inputs = ch_bam_bai
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .multiMap { meta, bam, bai, _meta_fasta, fasta, _meta_fai, fai ->
                bam: [meta, bam, bai]
                fasta: [meta, fasta, fai]
            }

        BAM_METHYLDACKEL(
            ch_methyldackel_inputs.bam,
            ch_methyldackel_inputs.fasta,
        )
        ch_bedgraph = BAM_METHYLDACKEL.out.methydackel_extract_bedgraph
        // channel: [ val(meta), [ bedgraph ]  ]
        ch_methylkit = BAM_METHYLDACKEL.out.methydackel_extract_methylkit
        // channel: [ val(meta), [ methylkit ] ]
        ch_mbias = BAM_METHYLDACKEL.out.methydackel_mbias
    }

    //
    // MODULE: Qualimap BamQC
    // skipped by default. to use run with `--run_qualimap` param.
    //
    if (params.run_qualimap) {
        QUALIMAP_BAMQC(
            ch_bam,
            params.bamqc_regions_file ? channel.fromPath(params.bamqc_regions_file, checkIfExists: true).toList() : [],
        )
        ch_qualimap = QUALIMAP_BAMQC.out.results
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
    }

    //
    // MODULE: Targeted sequencing analysis
    // skipped by default. to use run with `--run_targeted_sequencing` param.
    //
    if (params.run_targeted_sequencing) {
        if (params.taps || params.aligner == 'bwamem') {
            error("ERROR: --run_targeted_sequencing can't be running using rastair (methylation caller for TAPS) ")
        }
        if (!params.target_regions_file) {
            error("ERROR: --target_regions_file must be specified when using --run_targeted_sequencing")
        }
        TARGETED_SEQUENCING(
            ch_bedgraph,
            ch_coverage,
            channel.fromPath(params.target_regions_file, checkIfExists: true),
            ch_fasta,
            ch_fasta_index,
            ch_bam,
            ch_bai,
            ch_gzi,
            params.collecthsmetrics,
        )
        ch_versions = ch_versions.mix(TARGETED_SEQUENCING.out.versions)
    }

    //
    // MODULE: Preseq LCEXTRAP
    // skipped by default. to use run with `--run_preseq` param.
    //
    if (params.run_preseq) {
        PRESEQ_LCEXTRAP(
            ch_bam
        )
        ch_preseq = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'methylseq_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // Topic channel versions - collected below (after MULTIQC) and written to a
    // separate file, then merged into the main versions file on workflow completion.
    // MULTIQC's own version is mixed in explicitly because the nf-core module no
    // longer publishes to the `versions` topic (avoids a self-dependency hang).
    //
    ch_topic_versions = channel.topic("versions")

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_custom_methods_description = params.multiqc_methods_description
            ? file(params.multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true,
            )
        )

        if (params.run_qualimap) {
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect { it[1] }.ifEmpty([]))
        }
        if (params.run_preseq) {
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect { it[1] }.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc.ifEmpty([]))
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect { it[1] })
        }
        if (params.run_targeted_sequencing) {
            if (params.collecthsmetrics) {
                ch_multiqc_files = ch_multiqc_files.mix(TARGETED_SEQUENCING.out.picard_metrics.collect { it[1] }.ifEmpty([]))
            }
        }
        if (!params.skip_fastqc) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
        }

        // New nf-core MULTIQC (v4.0.2 template): single meta-based tuple input,
        // config passed as a list, logo/replace/sample resolved as values.
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'multiqc'],
                    files,
                    [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)] + (params.multiqc_config ? [file(params.multiqc_config, checkIfExists: true)] : []),
                    params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        // Keep the report path UNWRAPPED (flat) — PIPELINE_COMPLETION consumes it via getVal().
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> report }.toList()
        // MULTIQC now emits its version via `emit: versions` (topic-shaped tuple); fold it into the topic file.
        ch_topic_versions = ch_topic_versions.mix(MULTIQC.out.versions)
    }
    else {
        ch_multiqc_report = channel.empty()
    }

    //
    // Collate topic-channel versions (including MULTIQC) into a separate file.
    //
    ch_topic_versions
        .distinct()
        .filter { entry -> !(entry instanceof Path) }
        .map { process, tool, version ->
            def processName = process[process.lastIndexOf(':') + 1..-1]
            "${processName}:\n  ${tool}: ${version}"
        }
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_methylseq_topic_versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    bam            = ch_bam // channel: [ val(meta), path(bam) ]
    bai            = ch_bai // channel: [ val(meta), path(bai) ]
    rastair_mbias  = ch_rastair_mbias // channel: [ val(meta), path(rastair_mbias) ]
    rastair_call   = ch_rastair_call // channel: [ val(meta), path(rastair_call) ]
    methylkit      = ch_methylkit // channel: [ val(meta), path(methylkit) ]
    mbias          = ch_mbias // channel: [ val(meta), path(mbias) ]
    bedgraph       = ch_bedgraph // channel: [ val(meta), path(bedgraph) ]
    qualimap       = ch_qualimap // channel: [ val(meta), path(qualimap) ]
    preseq         = ch_preseq // channel: [ val(meta), path(preseq) ]
    multiqc_report = ch_multiqc_report // channel: [ path(multiqc_report.html )  ]
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
