//
// Alignment with BWA
//

nextflow.enable.types = true

include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BWA {
    take:
    ch_reads: Channel<Sample>
    val_index: Value<Path>
    sort_bam: Boolean
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>

    main:

    //
    // Map reads with BWA
    //
    ch_bam = BWA_MEM(ch_reads.combine(bwa_index: val_index, fasta: val_fasta), sort_bam)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_results = ch_bam
        .map { r -> record(meta: r.meta, align_bam: r.bam) }
        .join(BAM_SORT_STATS_SAMTOOLS(ch_bam, val_fasta, val_fasta_index), by: 'meta')

    emit:
    ch_results
}

record Sample {
    meta: Record
    reads: List<Path>
}
