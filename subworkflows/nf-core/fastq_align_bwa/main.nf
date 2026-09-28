//
// Alignment with BWA
//

nextflow.enable.types = true

include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

include { Sample } from '../../../utils/types.nf'

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
    ch_results = BAM_SORT_STATS_SAMTOOLS(ch_bam, val_fasta, val_fasta_index)

    emit:
    ch_results
}
