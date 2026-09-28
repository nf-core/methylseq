//
// Run SAMtools stats, flagstat and idxstats
//

nextflow.enable.types = true

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai: Channel<AlignedSample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>

    main:
    ch_stats = SAMTOOLS_STATS(ch_bam_bai.combine(fasta: val_fasta, fai: val_fasta_index))

    ch_flagstat = SAMTOOLS_FLAGSTAT(ch_bam_bai)

    ch_idxstats = SAMTOOLS_IDXSTATS(ch_bam_bai)

    ch_results = ch_stats
        .join(ch_flagstat, by: 'meta')
        .join(ch_idxstats, by: 'meta')

    emit:
    ch_results
}

record AlignedSample {
    meta: Record
    bam: Path
    bai: Path
}
