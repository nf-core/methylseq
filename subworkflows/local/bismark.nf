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
    genome // channel: [ val(meta), [ genome ] ]
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    /*
     * Generate bismark index if not supplied
     */

    // there might be indices from user input
    // branch them into different channels
    genome
        .branch{ meta, genome ->
            have_index: genome.containsKey('bismark_index')
                return [meta, genome.bismark_index]
            need_index: !genome.containsKey('bismark_index')
                return [meta, genome.fasta]
        }
        .set{ch_genome}

    // group by unique fastas, so that we only index each genome once
    ch_genome.need_index.groupTuple(by:1) | view | BISMARK_GENOMEPREPARATION

    // reverse groupTuple to restore the cardinality of the input channel
    // then join back with pre-existing indices
    bismark_index = BISMARK_GENOMEPREPARATION.out.index | transpose | mix(ch_genome.have_index)

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
