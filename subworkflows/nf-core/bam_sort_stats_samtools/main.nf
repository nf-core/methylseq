//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

nextflow.enable.types = true

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam: Channel<Alignment>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>

    main:
    ch_sorted = SAMTOOLS_SORT(ch_bam.combine(fasta: val_fasta, fai: val_fasta_index), '')
        .map { r -> record(meta: r.meta, bam: r.bam) }

    ch_bam_bai = ch_sorted.join(SAMTOOLS_INDEX(ch_sorted), by: 'meta')

    ch_stats = BAM_STATS_SAMTOOLS(ch_bam_bai, val_fasta, val_fasta_index)

    emit:
    ch_bam_bai.join(ch_stats, by: 'meta')
}

record Alignment {
    meta: Record
    bam: Path
}
