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

nextflow.preview.types = true

include { FASTA_INDEX_BISMARK_BWAMETH } from './subworkflows/nf-core/fasta_index_bismark_bwameth/main'
include { PIPELINE_INITIALISATION     } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { PIPELINE_COMPLETION         } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { getGenomeAttribute          } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { METHYLSEQ                   } from './workflows/methylseq/'

include { Sample                      } from './utils/types.nf'
include { MethylseqParams             } from './workflows/methylseq/'
include { MethylseqResult             } from './workflows/methylseq/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated file containing information about the samples in the experiment.
    input: String

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String

    /// MultiQC options

    // Custom config file to supply to MultiQC.
    multiqc_config: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: String?

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: String?

    /// Intermediate files

    // Save reference(s) to results directory
    save_reference: Boolean

    // Save aligned intermediates to results directory
    save_align_intermeds: Boolean

    // Bismark only - Save unmapped reads to FastQ files
    unmapped: Boolean

    // Save trimmed reads to results directory.
    save_trimmed: Boolean

    /// Reference options

    // Name of iGenomes reference.
    genome: String?

    // Path to FASTA genome file
    fasta: Path?

    // Path to Fasta index file.
    fasta_index: Path?

    // Path to a directory containing a Bismark reference index.
    bismark_index: Path?

    // bwameth index filename base
    bwameth_index: Path?

    /// Alignment options

    // Alignment tool to use.
    aligner: String = 'bismark'

    // Use BWA-MEM2 algorithm for BWA-Meth indexing and alignment.
    use_mem2: Boolean

    /// Library presets

    // Preset for working with PBAT libraries.
    pbat: Boolean

    // Turn on if dealing with MspI digested material.
    rrbs: Boolean

    // Run bismark in SLAM-seq mode.
    slamseq: Boolean

    // Preset for EM-seq libraries.
    em_seq: Boolean

    // Trimming preset for single-cell bisulfite libraries.
    single_cell: Boolean

    // Trimming preset for the Accel kit.
    accel: Boolean

    // Trimming preset for the Zymo kit.
    zymo: Boolean

    /// Trimming options

    // Trim bases from the 5' end of read 1 (or single-end reads).
    clip_r1: Integer = 0

    // Trim bases from the 5' end of read 2 (paired-end only).
    clip_r2: Integer = 0

    // Trim bases from the 3' end of read 1 AFTER adapter/quality trimming.
    three_prime_clip_r1: Integer = 0

    // Trim bases from the 3' end of read 2 AFTER adapter/quality trimming
    three_prime_clip_r2: Integer = 0

    // Trim bases below this quality value from the 3' end of the read, ignoring high-quality G bases
    nextseq_trim: Integer = 0

    // Discard reads that become shorter than INT because of either quality or adapter trimming.
    length_trim: Integer?

    // Skip presetting trimming parameters entirely
    skip_trimming_presets: Boolean

    /// Bismark options

    // Run alignment against all four possible strands.
    non_directional: Boolean

    // Output stranded cytosine report, following Bismark's bismark_methylation_extractor step.
    cytosine_report: Boolean

    // Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches).
    relax_mismatches: Boolean

    // 0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)
    num_mismatches: Float = 0.6

    // Specify a minimum read coverage to report a methylation call
    meth_cutoff: Integer?

    // Ignore read 2 methylation when it overlaps read 1
    no_overlap: Boolean = true

    // Ignore methylation in first n bases of 5' end of R1
    ignore_r1: Integer = 0

    // Ignore methylation in first n bases of 5' end of R2
    ignore_r2: Integer = 2

    // Ignore methylation in last n bases of 3' end of R1
    ignore_3prime_r1: Integer = 0

    // Ignore methylation in last n bases of 3' end of R2
    ignore_3prime_r2: Integer = 2

    // Supply a .gtf file containing known splice sites (bismark_hisat only).
    known_splices: String?

    // Allow soft-clipping of reads (potentially useful for single-cell experiments).
    local_alignment: Boolean

    // The minimum insert size for valid paired-end alignments.
    minins: Integer?

    // The maximum insert size for valid paired-end alignments.
    maxins: Integer?

    // Sample is NOMe-seq or NMT-seq. Runs coverage2cytosine.
    nomeseq: Boolean

    // Merges methylation calls for every strand into a single, context dependent file.
    comprehensive: Boolean

    /// bwa-meth options

    // Call methylation in all three CpG, CHG and CHH contexts.
    all_contexts: Boolean

    // Merges methylation metrics of the Cytosines in a given context.
    merge_context: Boolean

    // Specify a minimum read coverage for MethylDackel to report a methylation call.
    min_depth: Integer = 0

    // MethylDackel - ignore SAM flags
    ignore_flags: Boolean

    // Save files for use with methylKit
    methyl_kit: Boolean

    /// Qualimap options

    // A GFF or BED file containing the target regions which will be passed to Qualimap/Bamqc.
    bamqc_regions_file: String?

    /// Targeted sequencing options

    // A BED file containing the target regions
    target_regions_file: String?

    // Run Picard CollectHsMetrics in the targeted analysis
    collecthsmetrics: Boolean

    /// Skipping options

    // Skip read trimming.
    skip_trimming: Boolean

    // Skip deduplication step after alignment.
    skip_deduplication: Boolean

    // Skip FastQC
    skip_fastqc: Boolean

    // Skip MultiQC
    skip_multiqc: Boolean

    /// Run options

    // Run preseq/lcextrap tool
    run_preseq: Boolean

    // Run qualimap/bamqc tool
    run_qualimap: Boolean

    // Run advanced analysis for targeted methylation kits with enrichment of specific regions
    run_targeted_sequencing: Boolean

    // Email address for completion summary.
    email: String?

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Incoming hook URL for messaging service
    hook_url: String?

    // Display version and exit.
    version: Boolean

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true
}

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
    ch_samples: Channel<Sample>
    params_index: IndexParams
    params_methylseq: MethylseqParams

    main:

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    fasta = params_index.fasta ?: file(getGenomeAttribute('fasta', params))
    fasta_index = params_index.fasta_index ?: file(getGenomeAttribute('fasta_index', params))
    bismark_index = params_index.bismark_index ?: bismarkIndex(params_methylseq.aligner)
    bwameth_index = params_index.bwameth_index ?: bwamethIndex()

    indices = FASTA_INDEX_BISMARK_BWAMETH(
        fasta,
        fasta_index,
        bismark_index,
        bwameth_index,
        params_index.use_mem2,
        params_methylseq
    )

    //
    // WORKFLOW: Run pipeline
    //

    methylseq = METHYLSEQ (
        ch_samples,
        indices.fasta,
        indices.fasta_index,
        indices.bismark_index,
        indices.bwameth_index,
        params_methylseq
    )

    emit:
    fasta_index = indices.fasta_index
    bismark_index = indices.bismark_index
    bwameth_index = indices.bwameth_index
    results = methylseq.results
    bismark_summary = methylseq.bismark_summary
    reference_dict = methylseq.reference_dict
    intervallist = methylseq.intervallist
    multiqc_report = methylseq.multiqc_report

}

record IndexParams {
    fasta: Path?
    fasta_index: Path?
    bismark_index: Path?
    bwameth_index: Path?
    use_mem2: Boolean
}

def bismarkIndex(aligner: String) -> Path? {
    def indexPath = aligner == 'bismark_hisat'
        ? getGenomeAttribute('bismark_hisat2', params)
        : getGenomeAttribute('bismark', params)
    return indexPath ? file(indexPath) : null
}

def bwamethIndex() -> Path? {
    def indexPath = getGenomeAttribute('bwameth', params)
    return indexPath ? file(indexPath) : null
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
    ch_samples = PIPELINE_INITIALISATION (
        params.input,
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    methylseq = NFCORE_METHYLSEQ (
        ch_samples,
        params,
        params
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
        methylseq.multiqc_report
    )

    publish:
    fasta_index = methylseq.fasta_index
    bismark_index = methylseq.bismark_index
    bwameth_index = methylseq.bwameth_index
    samples = methylseq.results
    bismark_summary = methylseq.bismark_summary
    reference_dict = methylseq.reference_dict
    intervallist = methylseq.intervallist
    multiqc_report = methylseq.multiqc_report
}

output {
    fasta_index: Path {
        path "${params.aligner}/reference_genome"
        enabled params.save_reference
    }

    bismark_index: Path {
        path "${params.aligner}/reference_genome"
        enabled params.save_reference
    }

    bwameth_index: Path {
        path "${params.aligner}/reference_genome"
        enabled params.save_reference
    }

    samples: Channel<MethylseqResult> {
        path { r ->
            r.fastqc_html   >> "fastqc/"
            r.fastqc_zip    >> "fastqc/zips/"

            r.trim_reads    >> (params.save_trimmed ? "trimgalore/" : null)
            r.trim_log      >> "trimgalore/logs/"
            r.trim_unpaired >> (params.save_trimmed ? "trimgalore/" : null)
            r.trim_html     >> "trimgalore/fastqc/"
            r.trim_zip      >> "trimgalore/fastqc/zips/"

            r.bam   >> (params.save_align_intermeds ? "${params.aligner}/alignments/" : null)
            r.bai   >> (params.skip_deduplication ? "${params.aligner}/alignments/" : "${params.aligner}/deduplicated/")

            r.align_report                  >> "${params.aligner}/alignments/logs/"
            r.unmapped                      >> "${params.aligner}/alignments/unmapped/"
            r.dedup_report                  >> "${params.aligner}/deduplicated/logs/"
            r.coverage2cytosine_coverage    >> "bismark/coverage2cytosine/coverage/"
            r.coverage2cytosine_report      >> "bismark/coverage2cytosine/reports/"
            r.coverage2cytosine_summary     >> "bismark/coverage2cytosine/summaries/"
            r.methylation_bedgraph          >> "${params.aligner}/methylation_calls/bedGraph/"
            r.methylation_calls             >> "${params.aligner}/methylation_calls/methylation_calls/"
            r.methylation_coverage          >> "${params.aligner}/methylation_calls/methylation_coverage/"
            r.methylation_report            >> "${params.aligner}/methylation_calls/splitting_report/"
            r.methylation_mbias             >> "${params.aligner}/methylation_calls/mbias/"
            r.bismark_report                >> "${params.aligner}/reports/"

            r.samtools_flagstat                 >> "${params.aligner}/alignments/samtools_stats/"
            r.samtools_stats                    >> "${params.aligner}/alignments/samtools_stats/"
            r.methyldackel_extract_bedgraph     >> "methyldackel/"
            r.methyldackel_extract_methylkit    >> "methyldackel/"
            r.methyldackel_mbias                >> "methyldackel/mbias/"
            r.picard_metrics                    >> "${params.aligner}/deduplicated/picard_metrics/"

            r.qualimap_bamqc >> "${params.aligner}/qualimap/bamqc/"

            r.bedgraph_intersect >> (params.aligner == 'bismark' ? "bismark/methylation_calls/bedGraph/" : "methyldackel/")
            r.picard_hsmetrics >> "enrichment_metrics/"

            r.lc_extrap >> "${params.aligner}/preseq/"
            r.lc_log >> "${params.aligner}/preseq/log/"
        }
    }

    bismark_summary: Record {
        path "${params.aligner}/summary"
    }

    reference_dict: Record {
        path "${params.aligner}/reference_genome"
        enabled params.save_reference
    }

    intervallist: Record {
        path "enrichment_metrics"
    }

    multiqc_report: Path {
        path "multiqc/${params.aligner}"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
