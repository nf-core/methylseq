nextflow.enable.types = true

include { BISMARK_ALIGN                } from '../../../modules/nf-core/bismark/align/main'
include { BISMARK_DEDUPLICATE          } from '../../../modules/nf-core/bismark/deduplicate/main'
include { SAMTOOLS_SORT                } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX               } from '../../../modules/nf-core/samtools/index/main'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE    } from '../../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT               } from '../../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY              } from '../../../modules/nf-core/bismark/summary/main'

include { Sample } from '../../../utils/types.nf'

workflow FASTQ_ALIGN_DEDUP_BISMARK {
    take:
    ch_reads: Channel<Sample>
    val_fasta: Value<Path>
    val_bismark_index: Value<Path>
    skip_deduplication: Boolean // whether to deduplicate alignments
    cytosine_report: Boolean    // whether the run coverage2cytosine

    main:

    /*
     * Align with bismark
     */
    ch_alignment = BISMARK_ALIGN(
        ch_reads.combine(fasta: val_fasta, bismark_index: val_bismark_index)
    )
        .map { r -> r + record(align_bam: r.bam) }

    if (!skip_deduplication) {
        /*
        * Run deduplicate_bismark
        */
        ch_alignment_dedup = ch_alignment.join(BISMARK_DEDUPLICATE(ch_alignment), by: 'id')
    }
    else {
        ch_alignment_dedup = ch_alignment.map { r -> r + record(dedup_report: null) }
    }

    /*
     * MODULE: Run samtools sort on aligned or deduplicated bam
     */
    ch_bam = SAMTOOLS_SORT(ch_alignment_dedup.combine(fasta: val_fasta), '')

    /*
     * MODULE: Run samtools index on aligned or deduplicated bam
     */
    ch_bai = SAMTOOLS_INDEX(ch_bam)

    /*
     * Run bismark_methylation_extractor
     */
    ch_methylation = BISMARK_METHYLATIONEXTRACTOR(
        ch_alignment_dedup.combine(bismark_index: val_bismark_index)
    )

    /*
     * Run bismark coverage2cytosine
     */
    if (cytosine_report) {
        ch_coverage2cytosine = BISMARK_COVERAGE2CYTOSINE(
            ch_methylation.combine(fasta: val_fasta, bismark_index: val_bismark_index)
        )
    }
    else {
        ch_coverage2cytosine = channel.empty()
    }

    /*
     * Generate bismark sample reports
     */
    ch_bismark_report = BISMARK_REPORT(
        ch_alignment_dedup.join(ch_methylation, by: 'id')
    )

    /*
     * Collect per-sample results
     */
    ch_results = ch_alignment_dedup
        .join(ch_bam, by: 'id')
        .join(ch_bai, by: 'id')
        .join(ch_methylation, by: 'id')
        .join(ch_coverage2cytosine, by: 'id', remainder: true)
        .join(ch_bismark_report, by: 'id')

    /*
     * Generate bismark summary report
     */
    ch_bam_name = ch_alignment.map { r -> record(id: r.id, bam_name: r.bam.name) }

    val_bismark_summary_inputs = ch_alignment_dedup
        .join(ch_methylation, by: 'id')
        .join(ch_bam_name, by: 'id')
        .collect()
        .map { rs ->
            record(
                bam: rs.collect { r -> r.bam_name }.toSet(),
                align_report: rs.collect { r -> r.align_report }.toSet(),
                dedup_report: rs.collect { r -> r.dedup_report }.findAll { v -> v != null }.toSet(),
                methylation_report: rs.collect { r -> r.methylation_report }.toSet(),
                methylation_mbias: rs.collect { r -> r.methylation_mbias }.toSet()
            )
        }

    val_bismark_summary = BISMARK_SUMMARY(val_bismark_summary_inputs)

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = val_bismark_summary.flatMap { files -> files }
        .mix(ch_alignment_dedup.flatMap { r -> [r.align_report, r.dedup_report] }.filter { f -> f != null })
        .mix(ch_methylation.flatMap { r -> [r.methylation_report, r.methylation_mbias] })
        .mix(ch_bismark_report.flatMap { r -> r.bismark_report })

    emit:
    results         : Channel<BismarkResult> = ch_results
    bismark_summary : Value<Set<Path>>       = val_bismark_summary
    multiqc         : Channel<Path>          = ch_multiqc_files
}

record BismarkResult {
    id: String
    meta: Record
    bam: Path
    bai: Path
    align_bam: Path
    align_report: Path
    unmapped: Set<Path>
    dedup_report: Path?
    methylation_bedgraph: Path
    methylation_calls: Set<Path>
    methylation_coverage: Path
    methylation_report: Path
    methylation_mbias: Path
    coverage2cytosine_coverage: Path?
    coverage2cytosine_report: Path?
    coverage2cytosine_summary: Path?
    bismark_report: Set<Path>
}
