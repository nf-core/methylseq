nextflow.enable.types = true

include { BAM_SORT_STATS_SAMTOOLS       } from '../../nf-core/bam_sort_stats_samtools/main'
include { FASTQ_ALIGN_BWA               } from '../../nf-core/fastq_align_bwa/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES         } from '../../../modules/nf-core/picard/markduplicates/main'
include { PARABRICKS_FQ2BAM             } from '../../../modules/nf-core/parabricks/fq2bam/main'
include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_DEDUP_BWAMEM {
    take:
    ch_reads: Channel<Sample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>
    val_bwamem_index: Value<Path>
    skip_deduplication: Boolean // whether to deduplicate alignments
    use_gpu: Boolean            // whether to use GPU accelerated alignment
    output_fmt: String          // output format for parabricks fq2bam (e.g., 'bam' or 'cram')
    interval_file: List<Path>
    known_sites: List<Path>
    readgroups_args: String     // args for picard addorreplacereadgroups
    markduplicates_args: String // args for picard markduplicates

    main:
    /*
    Align with parabricks GPU enabled fq2bam implementation of bwa-mem
    */
    if (use_gpu) {
        ch_fq2bam = PARABRICKS_FQ2BAM(
            ch_reads.combine(fasta: val_fasta, bwa_index: val_bwamem_index),
            interval_file,
            known_sites,
            output_fmt,
        )
        ch_alignment = BAM_SORT_STATS_SAMTOOLS(ch_fq2bam, val_fasta, val_fasta_index)
    }
    else {
        ch_alignment = FASTQ_ALIGN_BWA(
            ch_reads,
            val_bwamem_index,
            true,
            val_fasta,
            val_fasta_index
        )
    }

    if (!skip_deduplication) {
        /*
         * Run Picard AddOrReplaceReadGroups to add read group (RG) to reads in bam file
         */
        ch_readgroups = PICARD_ADDORREPLACEREADGROUPS(
            ch_alignment
                .map { r -> record(meta: r.meta, bam: r.bam) }
                .combine(fasta: val_fasta, fai: val_fasta_index)
                .map { r -> r + record(args: readgroups_args) }
        )
        /*
         * Run Picard MarkDuplicates to mark duplicates
         */
        ch_markdup = PICARD_MARKDUPLICATES(
            ch_readgroups
                .combine(fasta: val_fasta, fai: val_fasta_index)
                .map { r -> r + record(args: markduplicates_args, prefix: "${r.meta.id}.markdup.sorted") }
        )
            .map { r -> record(meta: r.meta, bam: r.bam, picard_metrics: r.picard_metrics) }
        /*
         * Run samtools index on deduplicated alignment
         */
        ch_results = ch_alignment
            .join(ch_markdup, by: 'meta')
            .join(SAMTOOLS_INDEX(ch_markdup), by: 'meta')
    }
    else {
        ch_results = ch_alignment
            .map { r -> r + record(picard_metrics: null) }
    }

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_results.flatMap { r -> [r.picard_metrics, r.samtools_flagstat, r.samtools_stats, r.samtools_idxstats] }
        .filter { f -> f != null }

    emit:
    results : Channel<BwamemResult> = ch_results
    multiqc : Channel<Path>         = ch_multiqc_files
}

record BwamemResult {
    meta: Record
    bam: Path
    bai: Path
    align_bam: Path?
    samtools_flagstat: Path
    samtools_stats: Path
    samtools_idxstats: Path
    picard_metrics: Path?
}

record Sample {
    meta: Record
    reads: List<Path>
}
