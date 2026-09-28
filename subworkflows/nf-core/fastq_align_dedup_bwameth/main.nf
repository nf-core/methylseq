nextflow.enable.types = true

include { BWAMETH_ALIGN                                 } from '../../../modules/nf-core/bwameth/align/main'
include { PARABRICKS_FQ2BAMMETH                         } from '../../../modules/nf-core/parabricks/fq2bammeth/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS                                } from '../../../modules/nf-core/samtools/stats/main'
include { PICARD_MARKDUPLICATES                         } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'

include { Sample } from '../../../utils/types.nf'

workflow FASTQ_ALIGN_DEDUP_BWAMETH {
    take:
    ch_reads: Channel<Sample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>
    val_bwameth_index: Value<Path>
    skip_deduplication: Boolean // whether to deduplicate alignments
    use_gpu: Boolean            // whether to use GPU or CPU for bwameth alignment

    main:

    ch_align_inputs = ch_reads.combine(fasta: val_fasta, bwameth_index: val_bwameth_index)

    /*
     * Align with bwameth
     */
    if (use_gpu) {
        /*
        * Align with parabricks GPU enabled fq2bammeth implementation of bwameth
        */
        ch_alignment = PARABRICKS_FQ2BAMMETH(ch_align_inputs, [])
            .map { r -> record(id: r.id, meta: r.meta, bam: r.bam) }
    }
    else {
        /*
        * Align with CPU version of bwameth
        */
        ch_alignment = BWAMETH_ALIGN(ch_align_inputs)
    }

    /*
     * Sort raw output BAM
     */
    ch_sorted = SAMTOOLS_SORT(ch_alignment.combine(fasta: val_fasta, fai: val_fasta_index), '')
        .map { r -> record(id: r.id, meta: r.meta, bam: r.bam) }

    /*
     * Run samtools index on alignment
     */
    ch_sorted_bai = ch_sorted.join(SAMTOOLS_INDEX_ALIGNMENTS(ch_sorted), by: 'id')

    /*
     * Run samtools flagstat
     */
    ch_samtools_flagstat = SAMTOOLS_FLAGSTAT(ch_sorted_bai)

    /*
     * Run samtools stats
     */
    ch_samtools_stats = SAMTOOLS_STATS(ch_sorted_bai.combine(fasta: val_fasta, fai: val_fasta_index))

    if (!skip_deduplication) {
        /*
        * Run Picard MarkDuplicates
        */
        ch_markdup = PICARD_MARKDUPLICATES(ch_sorted.combine(fasta: val_fasta, fai: val_fasta_index))
            .map { r -> record(id: r.id, meta: r.meta, bam: r.bam, picard_metrics: r.picard_metrics) }

        /*
         * Run samtools index on deduplicated alignment
        */
        ch_alignment_final = ch_markdup.join(SAMTOOLS_INDEX_DEDUPLICATED(ch_markdup), by: 'id')
    }
    else {
        ch_alignment_final = ch_sorted_bai.map { r -> r + record(picard_metrics: null) }
    }

    ch_intermediates = ch_alignment
        .map { r -> record(id: r.id, align_bam: r.bam) }
        .join(ch_sorted_bai.map { r -> record(id: r.id, sorted_bam: r.bam, sorted_bai: r.bai) }, by: 'id')

    ch_results = ch_intermediates
        .join(ch_alignment_final, by: 'id')
        .join(ch_samtools_flagstat, by: 'id')
        .join(ch_samtools_stats, by: 'id')

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_results.flatMap { r -> [r.picard_metrics, r.samtools_flagstat, r.samtools_stats] }
        .filter { f -> f != null }

    emit:
    results : Channel<BwamethResult> = ch_results
    multiqc : Channel<Path>          = ch_multiqc_files
}

record BwamethResult {
    id: String
    meta: Record
    bam: Path
    bai: Path
    align_bam: Path
    sorted_bam: Path
    sorted_bai: Path
    samtools_flagstat: Path
    samtools_stats: Path
    picard_metrics: Path?
}
