/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.types = true

include { FASTQC                        } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE                    } from '../../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC                } from '../../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP               } from '../../modules/nf-core/preseq/lcextrap/main'
include { MULTIQC                       } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                     } from '../../modules/nf-core/cat/fastq/main'
include { WRITE_FILE                    } from '../../modules/local/writefile/main'
include { WRITE_FILE as WRITE_FILE_MULTIQC } from '../../modules/local/writefile/main'
include { FASTQ_ALIGN_DEDUP_BISMARK     } from '../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { FASTQ_ALIGN_DEDUP_BWAMETH     } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth/main'
include { FASTQ_ALIGN_DEDUP_BWAMEM      } from '../../subworkflows/nf-core/fastq_align_dedup_bwamem/main'
include { softwareVersionsToYAML        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { workflowSummaryMultiqc        } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { BAM_TAPS_CONVERSION           } from '../../subworkflows/nf-core/bam_taps_conversion'
include { BAM_METHYLDACKEL              } from '../../subworkflows/nf-core/bam_methyldackel/main'
include { TARGETED_SEQUENCING           } from '../../subworkflows/local/targeted_sequencing'
include { METHURATOR_GTESTIMATOR        } from '../../modules/nf-core/methurator/gtestimator/main'
include { METHURATOR_PLOT               } from '../../modules/nf-core/methurator/plot/main'

include { Sample                        } from '../../utils/types.nf'
include { toolArgs                      } from './args'
include { trimgaloreArgs                } from './args'
include { bismarkArgs                   } from './args'
include { bwamethArgs                   } from './args'
include { bwamemArgs                    } from './args'
include { methyldackelArgs              } from './args'
include { targetedSequencingArgs        } from './args'
include { qualimapArgs                  } from './args'
include { methuratorArgs                } from './args'
include { multiqcArgs                   } from './args'
include { runArgs                       } from './args'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {
    take:
    ch_samples: Channel<Sample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>?
    val_bismark_index: Value<Path>?
    val_bwameth_index: Value<Path>?
    val_bwamem_index: Value<Path>?
    params: MethylseqParams

    main:
    def ch_alignment: Channel<MethylseqResult> = channel.empty()
    def ch_aligner_mqc: Channel<Path> = channel.empty()
    def ch_bedgraph: Channel<SampleBedgraph> = channel.empty()
    def ch_methylation: Channel<SampleResult> = channel.empty()
    ch_multiqc_files = channel.empty()
    val_bismark_summary = null

    // sample metadata, joined back by id where needed
    ch_meta = ch_samples.map { s -> record(id: s.id, single_end: s.single_end, tool_args: s.tool_args) }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    ch_samples_multiple = ch_samples.filter { s -> s.reads.size() > (s.single_end ? 1 : 2) }
    ch_samples_single = ch_samples.filter { s -> s.reads.size() <= (s.single_end ? 1 : 2) }

    ch_fastq = ch_samples_multiple
        .join(CAT_FASTQ(ch_samples_multiple), by: 'id')
        .mix(ch_samples_single)

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        ch_fastqc = FASTQC(ch_fastq.map { s -> s + record(args: toolArgs('fastqc', s, params.args, '--quiet')) })
    }
    else {
        ch_fastqc = channel.empty()
    }

    //
    // MODULE: Run TrimGalore!
    //
    if (!params.skip_trimming) {
        ch_trimmed = TRIMGALORE(
            ch_fastq.map { s ->
                s + record(args: toolArgs('trimgalore', s, params.args, trimgaloreArgs(s, params)))
            }
        )
        ch_reads = ch_fastq
            .join(ch_trimmed, by: 'id')
            .map { r -> record(id: r.id, single_end: r.single_end, reads: r.trim_reads, tool_args: r.tool_args) }
    }
    else {
        ch_trimmed = channel.empty()
        ch_reads = ch_fastq
    }

    //
    // SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
    //

    if (params.taps && params.aligner != 'bwamem') {
        log.info("TAPS protocol detected and aligner is not 'bwamem'. We recommend using bwa-mem for TAPS protocol as it is optimized for this type of data.")
    }

    def use_gpu = workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1

    // Aligner: bismark or bismark_hisat
    if (params.aligner =~ /bismark/ && val_bismark_index) {
        //
        // Run Bismark alignment + downstream processing
        //
        bismark = FASTQ_ALIGN_DEDUP_BISMARK(
            ch_reads.map { r -> r + record(bismark_args: bismarkArgs(r, params)) },
            val_fasta,
            val_bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq,
        )
        ch_alignment = bismark.results
        ch_bedgraph = bismark.results.map { r -> record(id: r.id, bedgraphs: [r.methylation_bedgraph], coverage: r.methylation_coverage) }
        ch_aligner_mqc = bismark.multiqc
        val_bismark_summary = bismark.bismark_summary
    }
    else if (params.aligner == 'bwameth' && val_fasta_index && val_bwameth_index) {
        bwameth = FASTQ_ALIGN_DEDUP_BWAMETH(
            ch_reads.map { r -> r + record(bwameth_args: bwamethArgs(r, params.args)) },
            val_fasta,
            val_fasta_index,
            val_bwameth_index,
            params.skip_deduplication || params.rrbs,
            use_gpu,
        )
        ch_alignment = bwameth.results
        ch_aligner_mqc = bwameth.multiqc
    }
    else if (params.aligner == 'bwamem' && val_fasta_index && val_bwamem_index) {
        bwamem = FASTQ_ALIGN_DEDUP_BWAMEM(
            ch_reads.map { r -> r + record(bwamem_args: bwamemArgs(r, params.args)) },
            val_fasta,
            val_fasta_index,
            val_bwamem_index,
            params.skip_deduplication,
            use_gpu,
            'bam',
            [],
            [],
        )
        ch_alignment = bwamem.results
        ch_aligner_mqc = bwamem.multiqc
    }
    else {
        error("ERROR: Invalid aligner '${params.aligner}'. Valid options are: 'bismark', 'bismark_hisat', 'bwameth' or 'bwamem'.")
    }

    //
    // Subworkflow: Count positive mC->T conversion rates as a readout for DNA methylation
    //
    if ((params.taps || params.aligner == 'bwamem') && val_fasta_index) {
        ch_methylation = BAM_TAPS_CONVERSION(ch_alignment, val_fasta, val_fasta_index)
    }
    else if (!params.taps && params.aligner == 'bwameth' && val_fasta_index) {
        methyldackel = BAM_METHYLDACKEL(
            ch_alignment
                .join(ch_meta, by: 'id')
                .map { r -> r + record(methyldackel_args: methyldackelArgs(r, params)) },
            val_fasta,
            val_fasta_index
        )
        ch_methylation = methyldackel.results
        ch_bedgraph = methyldackel.results.map { r -> record(id: r.id, bedgraphs: r.methyldackel_bedgraph.toList(), coverage: null) }
    }

    //
    // MODULE: Qualimap BamQC
    // skipped by default. to use run with `--run_qualimap` param.
    //
    if (params.run_qualimap) {
        ch_qualimap = QUALIMAP_BAMQC(
            ch_alignment.join(ch_meta, by: 'id').map { r ->
                record(
                    id: r.id,
                    single_end: r.single_end,
                    bam: r.bam,
                    gff: params.bamqc_regions_file,
                    args: toolArgs('qualimap_bamqc', r, params.args, qualimapArgs(params))
                )
            }
        )
    }
    else {
        ch_qualimap = channel.empty()
    }

    //
    // MODULE: Targeted sequencing analysis
    // skipped by default. to use run with `--run_targeted_sequencing` param.
    //
    val_reference_dict = null
    val_intervallist = null
    if (params.run_targeted_sequencing && params.target_regions_file && val_fasta_index) {
        if (params.taps || params.aligner == 'bwamem') {
            error("ERROR: --run_targeted_sequencing can't be running using rastair (methylation caller for TAPS) ")
        }
        targeted_sequencing = TARGETED_SEQUENCING(
            ch_alignment
                .join(ch_bedgraph, by: 'id')
                .join(ch_meta, by: 'id')
                .map { r -> r + record(targeted_args: targetedSequencingArgs(r, params.args)) },
            channel.value(params.target_regions_file),
            val_fasta,
            val_fasta_index,
            params.collecthsmetrics,
            runArgs('picard_createsequencedictionary', params.args, ''),
            runArgs('picard_bedtointervallist', params.args, ''),
        )
        ch_targeted_sequencing = targeted_sequencing.results
        val_reference_dict = targeted_sequencing.reference_dict
        val_intervallist = targeted_sequencing.intervallist
    }
    else if (params.run_targeted_sequencing) {
        error("ERROR: --target_regions_file must be specified when using --run_targeted_sequencing")
    }
    else {
        ch_targeted_sequencing = channel.empty()
    }

    //
    // MODULE: Preseq LCEXTRAP
    // skipped by default. to use run with `--run_preseq` param.
    //
    if (params.run_preseq) {
        ch_preseq = PRESEQ_LCEXTRAP(
            ch_alignment.join(ch_meta, by: 'id').map { r ->
                record(id: r.id, single_end: r.single_end, bam: r.bam, args: toolArgs('preseq_lcextrap', r, params.args, ' -verbose -bam'))
            }
        )
    }
    else {
        ch_preseq = channel.empty()
    }

    //
    // MODULE: methurator gtestimator
    // skipped by default. to use run with `--run_methurator` param.
    //
    if (params.run_methurator && val_fasta_index) {
        if (params.taps || params.aligner == 'bwamem') {
            error("--run_methurator is not supported with the TAPS / bwa-mem workflow (methurator relies on MethylDackel).")
        }
        ch_methurator_gtestimator = METHURATOR_GTESTIMATOR(
            ch_alignment
                .combine(fasta: val_fasta, fai: val_fasta_index)
                .join(ch_meta, by: 'id')
                .map { r -> r + record(args: toolArgs('methurator_gtestimator', r, params.args, methuratorArgs(params))) }
        )
        ch_methurator = ch_methurator_gtestimator.join(METHURATOR_PLOT(ch_methurator_gtestimator), by: 'id')
    }
    else {
        ch_methurator = channel.empty()
    }

    //
    // Collect per-sample results
    //
    ch_results = ch_alignment
        .join(ch_fastqc, by: 'id', remainder: true)
        .join(ch_trimmed, by: 'id', remainder: true)
        .join(ch_methylation, by: 'id', remainder: true)
        .join(ch_qualimap, by: 'id', remainder: true)
        .join(ch_targeted_sequencing, by: 'id', remainder: true)
        .join(ch_preseq, by: 'id', remainder: true)
        .join(ch_methurator, by: 'id', remainder: true)

    //
    // Collate and save software versions
    //
    val_versions = softwareVersionsToYAML(channel.topic('versions'))
        .collect()
        .map { items ->
            record(name: 'nf_core_methylseq_software_mqc_versions.yml', items: items.toSorted(), newLine: true)
        }
    val_collated_versions = WRITE_FILE(val_versions)

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary = record(
            name: 'workflow_summary_mqc.yaml',
            items: [workflowSummaryMultiqc()]
        )

        multiqc_custom_methods_description = params.multiqc_methods_description
            ?: file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        methods_description = record(
            name: 'methods_description_mqc.yaml',
            items: [methodsDescriptionText(multiqc_custom_methods_description)]
        )

        ch_multiqc_files = ch_multiqc_files.mix(WRITE_FILE_MULTIQC(channel.of(workflow_summary, methods_description)))
        ch_multiqc_files = ch_multiqc_files.mix(val_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_qualimap.map { r -> r.qualimap_bamqc })
        ch_multiqc_files = ch_multiqc_files.mix(ch_preseq.map { r -> r.lc_log })
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc)
        ch_multiqc_files = ch_multiqc_files.mix(ch_trimmed.flatMap { r -> r.trim_log })
        ch_multiqc_files = ch_multiqc_files.mix(ch_targeted_sequencing.map { r -> r.picard_hsmetrics }.filter { f -> f != null })
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc.flatMap { r -> r.fastqc_zip })

        multiqc_default_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        multiqc_config = params.multiqc_config
            ? [multiqc_default_config, params.multiqc_config]
            : [multiqc_default_config]

        val_multiqc_inputs = ch_multiqc_files
            .collect()
            .map { files ->
                record(
                    multiqc_files: files.toSet(),
                    multiqc_config: multiqc_config,
                    multiqc_logo: params.multiqc_logo,
                    args: multiqcArgs(params)
                )
            }
        val_multiqc = MULTIQC(val_multiqc_inputs)
    }
    else {
        val_multiqc = null
    }

    emit:
    results         : Channel<MethylseqResult> = ch_results
    bismark_summary : Value<Set<Path>>?        = val_bismark_summary
    reference_dict  : Value<Path>?             = val_reference_dict
    intervallist    : Value<Path>?             = val_intervallist
    multiqc         : Value<MultiqcResult>?    = val_multiqc
    versions        : Channel<Path>            = val_collated_versions
}

record MethylseqParams {
    args: Map<String,String>
    slamseq: Boolean
    comprehensive: Boolean
    meth_cutoff: Integer?
    ignore_r1: Integer
    ignore_3prime_r1: Integer
    no_overlap: Boolean
    ignore_r2: Integer
    ignore_3prime_r2: Integer
    all_contexts: Boolean
    merge_context: Boolean
    ignore_flags: Boolean
    methyl_kit: Boolean
    min_depth: Integer
    genome: String?
    methurator_compute_ci: Boolean
    methurator_minimum_coverage: String?
    methurator_t_max: Integer?
    multiqc_title: String?
    aligner: String
    known_splices: Path?
    pbat: Boolean
    single_cell: Boolean
    non_directional: Boolean
    accel: Boolean
    zymo: Boolean
    em_seq: Boolean
    combined_index: Boolean
    local_alignment: Boolean
    unmapped: Boolean
    relax_mismatches: Boolean
    num_mismatches: Float
    minins: Integer?
    maxins: Integer?
    nextseq_trim: Integer
    length_trim: Integer?
    clip_r1: Integer
    clip_r2: Integer
    three_prime_clip_r1: Integer
    three_prime_clip_r2: Integer
    skip_trimming_presets: Boolean
    taps: Boolean
    skip_fastqc: Boolean
    skip_trimming: Boolean
    skip_deduplication: Boolean
    rrbs: Boolean
    cytosine_report: Boolean
    nomeseq: Boolean
    run_qualimap: Boolean
    bamqc_regions_file: Path?
    run_targeted_sequencing: Boolean
    target_regions_file: Path?
    collecthsmetrics: Boolean
    run_preseq: Boolean
    run_methurator: Boolean
    skip_multiqc: Boolean
    multiqc_config: Path?
    multiqc_logo: Path?
    multiqc_methods_description: Path?
}

record MultiqcResult {
    report: Path
    data: Path
    plots: Path?
}

record SampleResult {
    id: String
}

record SampleBedgraph {
    id: String
    bedgraphs: List<Path>
    coverage: Path?
}

record MethylseqResult {
    id: String
    single_end: Boolean

    // alignment (bismark / bwameth / bwamem)
    bam: Path
    bai: Path
    align_bam: Path?
    sorted_bam: Path?
    sorted_bai: Path?

    // fastqc
    fastqc_html: Set<Path>?
    fastqc_zip: Set<Path>?

    // trimgalore
    trim_reads: List<Path>?
    trim_log: List<Path>?
    trim_unpaired: List<Path>?
    trim_html: List<Path>?
    trim_zip: List<Path>?

    // bismark
    align_report: Path?
    unmapped: Set<Path>?
    dedup_report: Path?
    methylation_bedgraph: Path?
    methylation_calls: Set<Path>?
    methylation_coverage: Path?
    methylation_report: Path?
    methylation_mbias: Path?
    coverage2cytosine_coverage: Path?
    coverage2cytosine_report: Path?
    coverage2cytosine_summary: Path?
    bismark_report: Set<Path>?

    // bwameth / bwamem
    samtools_flagstat: Path?
    samtools_stats: Path?
    samtools_idxstats: Path?
    picard_metrics: Path?

    // methyldackel
    methyldackel_bedgraph: Set<Path>?
    methyldackel_methylkit: Set<Path>?
    methyldackel_mbias: Path?

    // rastair (taps)
    rastair_mbias: Path?
    rastair_mbias_pdf: Path?
    rastair_mbias_csv: Path?
    rastair_call: Path?
    rastair_methylkit: Path?

    // qualimap
    qualimap_bamqc: Path?

    // targeted sequencing
    bedgraph_intersect: Bag<Path>?
    coverage_intersect: Path?
    picard_hsmetrics: Path?

    // preseq
    lc_extrap: Path?
    lc_log: Path?

    // methurator
    methurator_summary: Path?
    methurator_plots: Set<Path>?
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
