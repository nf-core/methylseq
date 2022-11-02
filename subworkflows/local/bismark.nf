/*
 * bismark subworkflow
 */
include { BISMARK_GENOMEPREPARATION                   } from '../../modules/nf-core/bismark/genomepreparation/main'
include { BISMARK_ALIGN                               } from '../../modules/nf-core/bismark/align/main'
include { BISMARK_METHYLATIONEXTRACTOR                } from '../../modules/nf-core/bismark/methylationextractor/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED } from '../../modules/nf-core/samtools/sort/main'
include { BISMARK_DEDUPLICATE                         } from '../../modules/nf-core/bismark/deduplicate/main'
include { BISMARK_REPORT                              } from '../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                             } from '../../modules/nf-core/bismark/summary/main'

workflow BISMARK {
    take:
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    /*
     * Generate bismark index if not supplied
     */
    if (params.bismark_index) {
        bismark_index = file(params.bismark_index)
    } else {
        BISMARK_GENOMEPREPARATION(params.fasta)
        bismark_index = BISMARK_GENOMEPREPARATION.out.index
        ch_versions = ch_versions.mix(BISMARK_GENOMEPREPARATION.out.versions)
    }

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        reads,
        bismark_index
    )
    ch_versions = ch_versions.mix(BISMARK_ALIGN.out.versions)

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT_ALIGNED(
        BISMARK_ALIGN.out.bam,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_ALIGNED.out.versions)

    if (params.skip_deduplication || params.rrbs) {
        alignments = BISMARK_ALIGN.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    } else {
        /*
        * Run deduplicate_bismark
        */
        BISMARK_DEDUPLICATE( BISMARK_ALIGN.out.bam )

        alignments = BISMARK_DEDUPLICATE.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
        ch_versions = ch_versions.mix(BISMARK_DEDUPLICATE.out.versions)
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        bismark_index
    )
    ch_versions = ch_versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
        alignment_reports
            .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )
    ch_versions = ch_versions.mix(BISMARK_REPORT.out.versions)

    /*
     * Generate bismark summary report
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[2] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
    )
    ch_versions = ch_versions.mix(BISMARK_SUMMARY.out.versions)

    /*
     * MODULE: Run samtools sort
     */
    SAMTOOLS_SORT_DEDUPLICATED (
        alignments
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_DEDUPLICATED.out.versions)

    if (!params.skip_multiqc) {
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
    } else {
        multiqc_files = Channel.empty()
    }

    emit:
    bam        = SAMTOOLS_SORT_ALIGNED.out.bam        // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup      = SAMTOOLS_SORT_DEDUPLICATED.out.bam   // channel: [ val(meta), [ bam ] ] ## sorted, possibly deduplicated BAM
    mqc        = multiqc_files                        // path: *{html,txt}
    ch_versions                                       // path: *.version.txt
}
