/*
 * bismark subworkflow
 */

include { BISMARK_GENOMEPREPARATION    } from '../../modules/nf-core/modules/bismark/genomepreparation/main'
include { BISMARK_ALIGN                } from '../../modules/nf-core/modules/bismark/align/main'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../modules/nf-core/modules/bismark/methylationextractor/main'
include { SAMTOOLS_SORT                } from '../../modules/nf-core/modules/samtools/sort/main'
include { BISMARK_DEDUPLICATE          } from '../../modules/nf-core/modules/bismark/deduplicate/main'
include { BISMARK_REPORT               } from '../../modules/nf-core/modules/bismark/report/main'
include { BISMARK_SUMMARY              } from '../../modules/nf-core/modules/bismark/summary/main'

workflow BISMARK {
    take:
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    /*
     * Generate bismark index if not supplied
     */

    if (!params.genome || !params.genome.containsKey('bismarkIndex')) {
        BISMARK_GENOMEPREPARATION(params.fasta)
    }

    bismark_index = (!params.genome || !params.genome.containsKey('bismarkIndex')) ? BISMARK_GENOMEPREPARATION.out.index : params.genome.bismarkIndex

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        reads,
        bismark_index
    )

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
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        bismark_index
    )

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
    alignment_reports
        .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
        .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )

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

    /*
     * MODULE: Run samtools sort
     */
    SAMTOOLS_SORT (
        alignments
    )

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
            .set{ ch_multiqc_files }
    } else {
        ch_multiqc_files = Channel.empty()
    }

    /*
     * Collect modules Versions
     */
    SAMTOOLS_SORT.out.versions
        .mix(BISMARK_ALIGN.out.versions)
        .set{ versions }


    emit:
    bam              = BISMARK_ALIGN.out.bam          // channel: [ val(meta), [ bam ] ]
    dedup            = SAMTOOLS_SORT.out.bam          // channel: [ val(meta), [ bam ] ]

    mqc              = ch_multiqc_files               // path: *{html,txt}
    versions                                          // path: *.version.txt

}
