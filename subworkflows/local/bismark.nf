/*
 * bismark subworkflow
 */
include { BISMARK_ALIGN                               } from '../../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED } from '../../modules/nf-core/samtools/sort/main'
include { BISMARK_DEDUPLICATE                         } from '../../modules/nf-core/bismark/deduplicate/main'
include { BISMARK_METHYLATIONEXTRACTOR                } from '../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE                   } from '../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                              } from '../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                             } from '../../modules/nf-core/bismark/summary/main'

workflow BISMARK {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    bismark_index      // channel: /path/to/BismarkIndex/
    skip_deduplication // boolean: whether to deduplicate alignments
    cytosine_report    // boolean: whether the run coverage2cytosine

    main:
    versions = Channel.empty()


    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        reads,
        bismark_index
    )
    versions = versions.mix(BISMARK_ALIGN.out.versions)

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT_ALIGNED(
        BISMARK_ALIGN.out.bam,
    )
    versions = versions.mix(SAMTOOLS_SORT_ALIGNED.out.versions)

    if (skip_deduplication) {
        alignments = BISMARK_ALIGN.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    } else {
        /*
        * Run deduplicate_bismark
        */
        BISMARK_DEDUPLICATE( BISMARK_ALIGN.out.bam )

        alignments = BISMARK_DEDUPLICATE.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
        versions = versions.mix(BISMARK_DEDUPLICATE.out.versions)
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        bismark_index
    )
    versions = versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)


    /*
     * Run coverage2cytosine
     */
    if (cytosine_report) {
        BISMARK_COVERAGE2CYTOSINE (
            BISMARK_METHYLATIONEXTRACTOR.out.coverage,
            bismark_index
        )
        versions = versions.mix(BISMARK_COVERAGE2CYTOSINE.out.versions)
    }

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
        alignment_reports
            .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )
    versions = versions.mix(BISMARK_REPORT.out.versions)

    /*
     * Generate bismark summary report
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ it[1].name }.ifEmpty([]),
        alignment_reports.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[2] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
    )
    versions = versions.mix(BISMARK_SUMMARY.out.versions)

    /*
     * MODULE: Run samtools sort
     */
    SAMTOOLS_SORT_DEDUPLICATED (
        alignments
    )
    versions = versions.mix(SAMTOOLS_SORT_DEDUPLICATED.out.versions)

    /*
     * Collect MultiQC inputs
     */
    BISMARK_SUMMARY.out.summary.ifEmpty([])
        .mix(alignment_reports.collect{ it[1] })
        .mix(alignment_reports.collect{ it[2] })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] })
        .mix(BISMARK_REPORT.out.report.collect{ it[1] })
        .set{ multiqc_files }

    emit:
    bam        = SAMTOOLS_SORT_ALIGNED.out.bam        // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup      = SAMTOOLS_SORT_DEDUPLICATED.out.bam   // channel: [ val(meta), [ bam ] ] ## sorted, possibly deduplicated BAM
    mqc        = multiqc_files                        // path: *{html,txt}
    versions                                       // path: *.version.txt
}
