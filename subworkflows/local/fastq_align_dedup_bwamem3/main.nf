//
// Align TAPS reads with bwa-mem3 (`mem --meth=taps`), mark duplicates and collect alignment stats
//

include { BWAMEM3_MEM                            } from '../../../modules/nf-core/bwamem3/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGN } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_MARKDUPLICATES                  } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                     } from '../../nf-core/bam_stats_samtools/main'

workflow FASTQ_ALIGN_DEDUP_BWAMEM3 {
    take:
    ch_reads           // channel: [ val(meta), [ reads ] ]
    ch_fasta_fai       // channel: [ val(meta), fasta, fai ]
    ch_bwamem3_index   // channel: [ val(meta), [ bwa-mem3 --meth index ] ]
    skip_deduplication // boolean: whether to skip marking duplicates

    main:
    ch_picard_metrics = channel.empty()

    //
    // MODULE: Align and coordinate-sort with bwa-mem3
    //
    BWAMEM3_MEM(
        ch_reads,
        ch_bwamem3_index,
        ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        true,
    )
    SAMTOOLS_INDEX_ALIGN(BWAMEM3_MEM.out.aligned)
    ch_bam = BWAMEM3_MEM.out.aligned
    ch_bai = SAMTOOLS_INDEX_ALIGN.out.index

    //
    // MODULE: Mark duplicates
    //
    if (!skip_deduplication) {
        PICARD_MARKDUPLICATES(
            ch_bam,
            ch_fasta_fai,
        )
        SAMTOOLS_INDEX_DEDUP(PICARD_MARKDUPLICATES.out.bam)
        ch_bam = PICARD_MARKDUPLICATES.out.bam
        ch_bai = SAMTOOLS_INDEX_DEDUP.out.index
        ch_picard_metrics = PICARD_MARKDUPLICATES.out.metrics
    }

    //
    // SUBWORKFLOW: samtools stats, flagstat and idxstats on the final BAM
    //
    BAM_STATS_SAMTOOLS(
        ch_bam.join(ch_bai),
        ch_fasta_fai,
    )

    ch_multiqc_files = ch_picard_metrics
        .mix(BAM_STATS_SAMTOOLS.out.stats)
        .mix(BAM_STATS_SAMTOOLS.out.flagstat)
        .mix(BAM_STATS_SAMTOOLS.out.idxstats)
        .collect { _meta, file -> file }

    emit:
    bam     = ch_bam           // channel: [ val(meta), bam ]
    bai     = ch_bai           // channel: [ val(meta), bai ]
    multiqc = ch_multiqc_files // channel: [ files ]
}
