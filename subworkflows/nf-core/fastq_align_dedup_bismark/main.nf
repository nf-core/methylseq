include { BISMARK_ALIGN as BISMARK_ALIGN_HOST } from '../../../modules/nf-core/bismark/align/main'
include { BISMARK_ALIGN as BISMARK_ALIGN_SPIKE } from '../../../modules/nf-core/bismark/align/main'

include { BISMARK_DEDUPLICATE } from '../../../modules/nf-core/bismark/deduplicate/main'

include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

include { BISMARK_METHYLATIONEXTRACTOR as BISMARK_METHYLATIONEXTRACTOR_HOST } \
    from '../../../modules/nf-core/bismark/methylationextractor/main'

include { BISMARK_METHYLATIONEXTRACTOR as BISMARK_METHYLATIONEXTRACTOR_SPIKE } \
    from '../../../modules/nf-core/bismark/methylationextractor/main'

include { BISMARK_COVERAGE2CYTOSINE } \
    from '../../../modules/nf-core/bismark/coverage2cytosine/main'

include { BISMARK_REPORT } \
    from '../../../modules/nf-core/bismark/report/main'

include { BISMARK_SUMMARY } \
    from '../../../modules/nf-core/bismark/summary/main'


workflow FASTQ_ALIGN_DEDUP_BISMARK {

    take:
    ch_reads
    ch_fasta
    ch_bismark_index

    ch_spikein_fasta
    ch_spikein_bismark_index

    skip_deduplication
    cytosine_report
    spikein

    main:

    ch_alignments                 = channel.empty()
    ch_alignment_reports          = channel.empty()

    ch_methylation_bedgraph       = channel.empty()
    ch_methylation_calls          = channel.empty()
    ch_methylation_coverage       = channel.empty()
    ch_methylation_report         = channel.empty()
    ch_methylation_mbias          = channel.empty()

    ch_coverage2cytosine_coverage = channel.empty()
    ch_coverage2cytosine_report   = channel.empty()
    ch_coverage2cytosine_summary  = channel.empty()

    ch_bismark_report             = channel.empty()
    ch_bismark_summary            = channel.empty()

    ch_multiqc_files              = channel.empty()
    ch_versions                   = channel.empty()

    // Spike-in channels
    ch_spike_bam                  = channel.empty()
    ch_spike_alignment_reports    = channel.empty()
    ch_spike_methylation_report   = channel.empty()
    ch_spike_mbias                = channel.empty()
    ch_spike_coverage             = channel.empty()


    /*
     * HOST ALIGNMENT
     */
    BISMARK_ALIGN_HOST (
        ch_reads,
        ch_fasta,
        ch_bismark_index
    )

    ch_alignments =
        BISMARK_ALIGN_HOST.out.bam

    ch_alignment_reports =
        BISMARK_ALIGN_HOST.out.report.map { meta, report ->
            [ meta, report, [] ]
        }


    /*
     * SPIKE-IN ALIGNMENT
     */
    if (spikein) {

        BISMARK_ALIGN_SPIKE (
            ch_reads,
            ch_spikein_fasta,
            ch_spikein_bismark_index
        )

        ch_spike_bam =
            BISMARK_ALIGN_SPIKE.out.bam

        ch_spike_alignment_reports =
            BISMARK_ALIGN_SPIKE.out.report

        ch_versions =
            ch_versions.mix(
                BISMARK_ALIGN_SPIKE.out.versions
            )
    }


    /*
     * HOST DEDUPLICATION
     */
    if (!skip_deduplication) {

        BISMARK_DEDUPLICATE (
            BISMARK_ALIGN_HOST.out.bam
        )

        ch_alignments =
            BISMARK_DEDUPLICATE.out.bam

        ch_alignment_reports =
            BISMARK_ALIGN_HOST.out.report.join(
                BISMARK_DEDUPLICATE.out.report
            )

        ch_versions =
            ch_versions.mix(
                BISMARK_DEDUPLICATE.out.versions
            )
    }


    /*
     * SORT BAM
     */
    SAMTOOLS_SORT (
        ch_alignments,
        [[:],[]],
        ''
    )


    /*
     * INDEX BAM
     */
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )

    ch_versions =
        ch_versions.mix(
            SAMTOOLS_INDEX.out.versions
        )


    /*
     * HOST METHYLATION EXTRACTION
     */
    BISMARK_METHYLATIONEXTRACTOR_HOST (
        ch_alignments,
        ch_bismark_index
    )

    ch_methylation_bedgraph =
        BISMARK_METHYLATIONEXTRACTOR_HOST.out.bedgraph

    ch_methylation_calls =
        BISMARK_METHYLATIONEXTRACTOR_HOST.out.methylation_calls

    ch_methylation_coverage =
        BISMARK_METHYLATIONEXTRACTOR_HOST.out.coverage

    ch_methylation_report =
        BISMARK_METHYLATIONEXTRACTOR_HOST.out.report

    ch_methylation_mbias =
        BISMARK_METHYLATIONEXTRACTOR_HOST.out.mbias

    ch_versions =
        ch_versions.mix(
            BISMARK_METHYLATIONEXTRACTOR_HOST.out.versions
        )


    /*
     * SPIKE-IN METHYLATION EXTRACTION
     */
    if (spikein) {

        BISMARK_METHYLATIONEXTRACTOR_SPIKE (
            ch_spike_bam,
            ch_spikein_bismark_index
        )

        ch_spike_methylation_report =
            BISMARK_METHYLATIONEXTRACTOR_SPIKE.out.report

        ch_spike_mbias =
            BISMARK_METHYLATIONEXTRACTOR_SPIKE.out.mbias

        ch_spike_coverage =
            BISMARK_METHYLATIONEXTRACTOR_SPIKE.out.coverage

        ch_versions =
            ch_versions.mix(
                BISMARK_METHYLATIONEXTRACTOR_SPIKE.out.versions
            )
    }


    /*
     * COVERAGE2CYTOSINE
     */
    if (cytosine_report) {

        BISMARK_COVERAGE2CYTOSINE (
            ch_methylation_coverage,
            ch_fasta,
            ch_bismark_index
        )

        ch_coverage2cytosine_coverage =
            BISMARK_COVERAGE2CYTOSINE.out.coverage

        ch_coverage2cytosine_report =
            BISMARK_COVERAGE2CYTOSINE.out.report

        ch_coverage2cytosine_summary =
            BISMARK_COVERAGE2CYTOSINE.out.summary

        ch_versions =
            ch_versions.mix(
                BISMARK_COVERAGE2CYTOSINE.out.versions
            )
    }


    /*
     * HOST REPORTS
     */
    BISMARK_REPORT (
        ch_alignment_reports
            .join(ch_methylation_report)
            .join(ch_methylation_mbias)
    )

    ch_bismark_report =
        BISMARK_REPORT.out.report

    ch_versions =
        ch_versions.mix(
            BISMARK_REPORT.out.versions
        )


    /*
     * SUMMARY
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN_HOST.out.bam.collect { _meta, bam -> bam.name },
        ch_alignment_reports.collect { _meta, align_report, _dedup_report -> align_report },
        ch_alignment_reports.collect { _meta, _align_report, dedup_report -> dedup_report }.ifEmpty([]),
        ch_methylation_report.collect { _meta, report -> report },
        ch_methylation_mbias.collect { _meta, mbias -> mbias }
    )

    ch_bismark_summary =
        BISMARK_SUMMARY.out.summary

    ch_versions =
        ch_versions.mix(
            BISMARK_SUMMARY.out.versions
        )


    /*
     * MULTIQC
     */
    ch_multiqc_files =
        ch_bismark_summary
            .mix(ch_alignment_reports.collect { _meta, align_report, _dedup_report -> align_report })
            .mix(ch_alignment_reports.collect { _meta, _align_report, dedup_report -> dedup_report })
            .mix(ch_methylation_report.collect { _meta, report -> report })
            .mix(ch_methylation_mbias.collect { _meta, mbias -> mbias })
            .mix(ch_bismark_report.collect { _meta, report -> report })

    if (spikein) {
        ch_multiqc_files =
            ch_multiqc_files
                .mix(ch_spike_methylation_report.collect { _meta, report -> report })
                .mix(ch_spike_mbias.collect { _meta, mbias -> mbias })
    }

    emit:

    bam                        = SAMTOOLS_SORT.out.bam
    bai                        = SAMTOOLS_INDEX.out.bai

    coverage2cytosine_coverage = ch_coverage2cytosine_coverage
    coverage2cytosine_report   = ch_coverage2cytosine_report
    coverage2cytosine_summary  = ch_coverage2cytosine_summary

    methylation_bedgraph       = ch_methylation_bedgraph
    methylation_calls          = ch_methylation_calls
    methylation_coverage       = ch_methylation_coverage
    methylation_report         = ch_methylation_report
    methylation_mbias          = ch_methylation_mbias

    spike_bam                  = ch_spike_bam
    spike_coverage             = ch_spike_coverage
    spike_methylation_report   = ch_spike_methylation_report
    spike_mbias                = ch_spike_mbias

    bismark_report             = ch_bismark_report
    bismark_summary            = ch_bismark_summary

    multiqc                    = ch_multiqc_files
    versions                   = ch_versions
}

