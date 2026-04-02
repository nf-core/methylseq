
nextflow.preview.types = true

include { BWAMETH_ALIGN                                 } from '../../../modules/nf-core/bwameth/align/main'
include { PARABRICKS_FQ2BAMMETH                         } from '../../../modules/nf-core/parabricks/fq2bammeth/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS                                } from '../../../modules/nf-core/samtools/stats/main'
include { PICARD_MARKDUPLICATES                         } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'
include { METHYLDACKEL_EXTRACT                          } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS                            } from '../../../modules/nf-core/methyldackel/mbias/main'

include { Sample } from '../../../utils/types.nf'

workflow FASTQ_ALIGN_DEDUP_BWAMETH {

    take:
    ch_reads: Channel<Sample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>
    val_bwameth_index: Value<Path>
    skip_deduplication: Boolean
    use_gpu: Boolean

    main:

    /*
     * Align with bwameth
     */
    if (use_gpu) {
        /*
        * Align with parabricks GPU enabled fq2bammeth implementation of bwameth
        */
        ch_alignment = PARABRICKS_FQ2BAMMETH (
            ch_reads.combine(fasta: val_fasta, bwameth_index: val_bwameth_index)
        )
    } else {
        /*
        * Align with CPU version of bwameth
        */
        ch_alignment = BWAMETH_ALIGN (
            ch_reads.combine(fasta: val_fasta, bwameth_index: val_bwameth_index)
        )
    }

    /*
     * Sort raw output BAM
     */
    ch_alignment = SAMTOOLS_SORT( ch_alignment )

    /*
     * Run samtools index on alignment
     */
    ch_alignment_index = SAMTOOLS_INDEX_ALIGNMENTS(
        ch_alignment.map { r -> record(id: r.id, input: r.bam) }
    )
    ch_alignment = ch_alignment.join(ch_alignment_index, by: 'id')

    /*
     * Run samtools flagstat
     */
    ch_samtools_flagstat = SAMTOOLS_FLAGSTAT( ch_alignment )

    /*
     * Run samtools stats
     */
    ch_samtools_stats = SAMTOOLS_STATS(
        ch_alignment.map { r -> record(id: r.id, input: r.bam, input_index: r.bai) }
    )

    if (!skip_deduplication) {
        /*
        * Run Picard MarkDuplicates
        */
        ch_picard = PICARD_MARKDUPLICATES(
            ch_alignment.combine(fasta: val_fasta, fasta_index: val_fasta_index)
        )
        /*
         * Run samtools index on deduplicated alignment
        */
        ch_alignment_index_dedup = SAMTOOLS_INDEX_DEDUPLICATED(
            ch_picard.map { r -> record(id: r.id, input: r.bam) }
        )
        ch_alignment = ch_alignment
            .join(ch_picard, by: 'id')
            .join(ch_alignment_index_dedup, by: 'id')
    }

    /*
     * Extract per-base methylation and plot methylation bias
     */

    ch_methydackel_extract = METHYLDACKEL_EXTRACT (
        ch_alignment.combine(fasta: val_fasta, fasta_index: val_fasta_index)
    )

    ch_methydackel_mbias = METHYLDACKEL_MBIAS (
        ch_alignment.combine(fasta: val_fasta, fasta_index: val_fasta_index)
    )

    ch_results = ch_alignment
        .join(ch_samtools_flagstat, by: 'id')
        .join(ch_samtools_stats, by: 'id')
        .join(ch_methydackel_extract, by: 'id')
        .join(ch_methydackel_mbias, by: 'id')
        .join(ch_picard, by: 'id', remainder: true)

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = channel.empty()
        .mix( ch_picard.map { r -> r.picard_metrics } )
        .mix( ch_samtools_flagstat.map { r -> r.samtools_flagstat } )
        .mix( ch_samtools_stats.map { r -> r.samtools_stats } )
        .mix( ch_methydackel_extract.map { r -> r.methydackel_bedgraph } )
        .mix( ch_methydackel_mbias.map { r -> r.methyldackel_mbias } )

    emit:
    results: Channel<BwamethResult> = ch_results
    multiqc: Channel<Path> = ch_multiqc_files
}


record BwamethResult {
    id: String
    single_end: Boolean
    bam: Path
    bai: Path
    samtools_flagstat: Path
    samtools_stats: Path
    methyldackel_extract_bedgraph: Path
    methyldackel_extract_methylkit: Path
    methyldackel_mbias: Path
    picard_metrics: Path
}
