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

nextflow.enable.types = true

include { FASTA_INDEX_METHYLSEQ     } from './subworkflows/nf-core/fasta_index_methylseq/main'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { getGenomeAttribute        } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { METHYLSEQ                 } from './workflows/methylseq/'
include { MethylseqParams           } from './workflows/methylseq/'
include { Sample                    } from './workflows/methylseq/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    /// Input/output options

    // Path to comma-separated file containing information about the samples in the experiment.
    input: Path?

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String?

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    /// Save intermediate files

    // Save reference(s) to results directory
    save_reference: Boolean

    // Save aligned intermediates to results directory
    save_align_intermeds: Boolean

    // Bismark only - Save unmapped reads to FastQ files
    unmapped: Boolean

    // Save trimmed reads to results directory.
    save_trimmed: Boolean

    /// Reference genome options

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

    // Path to the BWA-MEM index filename base
    bwamem_index: Path?

    /// Alignment options

    // Alignment tool to use.
    aligner: String = 'bismark'

    // Use BWA-MEM2 algorithm for BWA-Meth indexing and alignment.
    use_mem2: Boolean

    /// Special library types

    // Preset for working with TET-assisted pyridine borane sequencing (TAPS) libraries.
    taps: Boolean

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

    /// Adapter Trimming

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
    known_splices: Path?

    // Allow soft-clipping of reads (potentially useful for single-cell experiments).
    local_alignment: Boolean

    // Align against a single combined (CT+GA) Bismark index instead of the classic per-strand model.
    combined_index: Boolean

    // The minimum insert size for valid paired-end alignments.
    minins: Integer?

    // The maximum insert size for valid paired-end alignments.
    maxins: Integer?

    // Sample is NOMe-seq or NMT-seq. Runs coverage2cytosine.
    nomeseq: Boolean

    // Merges methylation calls for every strand into a single, context dependent file.
    comprehensive: Boolean

    /// MethylDackel options

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

    /// rastair options

    // Nucleotides to exclude from methylation calling on Original Top (strand) reads, following this pattern: R1 start, R1 end, R2 start, R2 end
    trim_OT: String = '0,0,10,0'

    // Nucleotides to exclude from methylation calling on Original Bottom (strand) reads, following this pattern: R1 start, R1 end, R2 start, R2 end
    trim_OB: String = '0,0,10,0'

    /// methurator options

    // Minimum CpGs coverage to consider for the saturation analysis. Can be a single integer or a list (e.g. 1,3,5).
    methurator_minimum_coverage: String = '1,10,15'

    // Compute confidence intervals using bootstrap replicates.
    methurator_compute_ci: Boolean

    // Maximum extrapolation factor.
    methurator_t_max: Integer = 10

    /// Qualimap Options

    // A GFF or BED file containing the target regions which will be passed to Qualimap/Bamqc.
    bamqc_regions_file: Path?

    /// Targeted Sequencing Analysis Options

    // A BED file containing the target regions
    target_regions_file: Path?

    // Run Picard CollectHsMetrics in the targeted analysis
    collecthsmetrics: Boolean

    /// Skip pipeline steps

    // Skip read trimming.
    skip_trimming: Boolean

    // Skip deduplication step after alignment.
    skip_deduplication: Boolean

    // Skip FastQC
    skip_fastqc: Boolean

    // Skip MultiQC
    skip_multiqc: Boolean

    /// Run pipeline steps

    // Run preseq/lcextrap tool
    run_preseq: Boolean

    // Run methurator/gtestimator tool
    run_methurator: Boolean

    // Run qualimap/bamqc tool
    run_qualimap: Boolean

    // Run advanced analysis for targeted methylation kits with enrichment of specific regions
    run_targeted_sequencing: Boolean

    /// Generic options

    // Display version and exit.
    version: Boolean

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Display the help message.
    help: String?

    // Display the full detailed help message.
    help_full: Boolean

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean
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
    // Initialize reference files from params or iGenomes
    //
    fasta         = params_index.fasta         ?: genomeFile('fasta')
    fasta_index   = params_index.fasta_index   ?: genomeFile('fasta_index')
    bismark_index = params_index.bismark_index ?: genomeFile(params_methylseq.aligner == 'bismark_hisat' ? 'bismark_hisat2' : 'bismark')
    bwameth_index = params_index.bwameth_index ?: genomeFile('bwameth')
    bwamem_index  = params_index.bwamem_index  ?: genomeFile('bwa')

    if (!fasta) {
        error("ERROR: A reference genome must be provided with --fasta or --genome")
    }

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    indices = FASTA_INDEX_METHYLSEQ(
        fasta,
        fasta_index,
        bismark_index,
        bwameth_index,
        bwamem_index,
        params_methylseq.aligner,
        params_methylseq.collecthsmetrics,
        params_methylseq.run_methurator,
        params_index.use_mem2
    )

    //
    // WORKFLOW: Run pipeline
    //
    methylseq = METHYLSEQ(
        ch_samples,
        indices.fasta,
        indices.fasta_index,
        indices.bismark_index,
        indices.bwameth_index,
        indices.bwamem_index,
        params_methylseq
    )

    emit:
    methylseq.multiqc_report
}

record IndexParams {
    fasta: Path?
    fasta_index: Path?
    bismark_index: Path?
    bwameth_index: Path?
    bwamem_index: Path?
    use_mem2: Boolean
}

def genomeFile(attribute: String) -> Path? {
    def path = getGenomeAttribute(attribute)
    return path ? file(path) : null
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
    ch_samples = PIPELINE_INITIALISATION(
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
    val_multiqc_report = NFCORE_METHYLSEQ(
        ch_samples,
        params,
        params
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        val_multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
