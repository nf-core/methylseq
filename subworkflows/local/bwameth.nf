/*
 * bwameth subworkflow
 */
include { SAMTOOLS_STATS                                } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX                                } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES                         } from '../../modules/nf-core/picard/markduplicates/main'
include { BWAMETH_INDEX                                 } from '../../modules/nf-core/bwameth/index/main'
include { BWAMETH_ALIGN                                 } from '../../modules/nf-core/bwameth/align/main'
include { METHYLDACKEL_EXTRACT                          } from '../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS                            } from '../../modules/nf-core/methyldackel/mbias/main'

workflow BWAMETH {
    take:
    reads  // channel: [ val(meta), [ reads ] ]

    main:

    versions = Channel.empty()

    /*
     * Generate bwameth index if not supplied
     */
    if (params.bwa_meth_index) {
        bwameth_index = file(params.bwa_meth_index)
    } else {
        BWAMETH_INDEX(params.fasta)
        bwameth_index = BWAMETH_INDEX.out.index
        versions = versions.mix(BWAMETH_INDEX.out.versions)
    }

    /*
     * Generate fasta index if not supplied
     */
    if (params.fasta_index) {
        fasta_index = file(params.fasta_index)
    } else {
        SAMTOOLS_FAIDX([[:], params.fasta])
        fasta_index = SAMTOOLS_FAIDX.out.fai.map{ return(it[1])}
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    /*
     * Align with bwameth
     */
    BWAMETH_ALIGN (
        reads,
        bwameth_index
    )
    versions = versions.mix(BWAMETH_ALIGN.out.versions)

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT (BWAMETH_ALIGN.out.bam)
    versions = versions.mix(SAMTOOLS_SORT.out.versions)

    /*
     * Run samtools index on alignment
     */
    SAMTOOLS_INDEX_ALIGNMENTS (SAMTOOLS_SORT.out.bam)
    versions = versions.mix(SAMTOOLS_INDEX_ALIGNMENTS.out.versions)

    /*
     * Run samtools flagstat and samtools stats
     */
    SAMTOOLS_FLAGSTAT( BWAMETH_ALIGN.out.bam.join(SAMTOOLS_INDEX_ALIGNMENTS.out.bai) )
    SAMTOOLS_STATS( BWAMETH_ALIGN.out.bam.join(SAMTOOLS_INDEX_ALIGNMENTS.out.bai), [] )
    versions = versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    versions = versions.mix(SAMTOOLS_STATS.out.versions)

    if (params.skip_deduplication || params.rrbs) {
        alignments = SAMTOOLS_SORT.out.bam
        bam_index = SAMTOOLS_INDEX_ALIGNMENTS.out.bai
        picard_metrics = Channel.empty()
        picard_version = Channel.empty()
    } else {
        /*
        * Run Picard MarkDuplicates
        */
        PICARD_MARKDUPLICATES (
            SAMTOOLS_SORT.out.bam,
            params.fasta,
            fasta_index
        )
        /*
         * Run samtools index on deduplicated alignment
        */
        SAMTOOLS_INDEX_DEDUPLICATED (PICARD_MARKDUPLICATES.out.bam)

        alignments = PICARD_MARKDUPLICATES.out.bam
        bam_index = SAMTOOLS_INDEX_DEDUPLICATED.out.bai
        picard_metrics = PICARD_MARKDUPLICATES.out.metrics
        picard_version = PICARD_MARKDUPLICATES.out.versions
        versions = versions.mix(PICARD_MARKDUPLICATES.out.versions)
    }


    /*
     * Extract per-base methylation and plot methylation bias
     */

    METHYLDACKEL_EXTRACT(
        alignments.join(bam_index),
        params.fasta,
        fasta_index
    )
    METHYLDACKEL_MBIAS(
        alignments.join(bam_index),
        params.fasta,
        fasta_index
    )
    versions = versions.mix(METHYLDACKEL_EXTRACT.out.versions)
    versions = versions.mix(METHYLDACKEL_MBIAS.out.versions)

    /*
     * Collect MultiQC inputs
     */
    if (!params.skip_multiqc) {
        picard_metrics.collect{ it[1] }
            .mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{ it[1] })
            .mix(SAMTOOLS_STATS.out.stats.collect{ it[1] })
            .mix(METHYLDACKEL_EXTRACT.out.bedgraph.collect{ it[1] })
            .mix(METHYLDACKEL_MBIAS.out.txt.collect{ it[1] })
            .set{ multiqc_files }
    } else {
        multiqc_files = Channel.empty()
    }

    emit:
    bam                  = SAMTOOLS_SORT.out.bam  // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup                = alignments             // channel: [ val(meta), [ bam ] ]  ## sorted, possibly deduplicated BAM
    mqc                  = multiqc_files          // path: *{html,txt}
    versions                                   // path: *.version.txt
}
