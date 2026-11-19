//
// UMI libraries with bwa-mem3: fgumi runall extracts the UMIs, aligns, restores the UMI tags and sorts by template
// coordinate in one process; fgumi dedup then marks duplicates per UMI family and strand
//

include { SAMTOOLS_DICT      } from '../../../modules/nf-core/samtools/dict/main'
include { FGUMI_RUNALL       } from '../../../modules/local/fgumi_runall/main'
include { FGUMI_DEDUP        } from '../../../modules/nf-core/fgumi/dedup/main'
include { FGUMI_SORT         } from '../../../modules/nf-core/fgumi/sort/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow FASTQ_ALIGN_DEDUP_UMI_BWAMEM3 {
    take:
    ch_reads           // channel: [ val(meta), [ reads ] ]
    ch_fasta_fai       // channel: [ val(meta), fasta, fai ]
    ch_bwamem3_index   // channel: [ val(meta), [ bwa-mem3 index ] ]
    skip_deduplication // boolean: whether to skip marking duplicates

    main:
    ch_dedup_metrics = channel.empty()

    // One reference for every sample
    ch_reference = ch_fasta_fai.first()
    SAMTOOLS_DICT(ch_reference.map { meta, fasta, _fai -> [meta, fasta] })

    FGUMI_RUNALL(
        ch_reads.map { meta, reads -> [meta, reads, meta.id] },
        ch_reference.join(SAMTOOLS_DICT.out.dict).first(),
        ch_bwamem3_index,
        'extract',
        'sort',
    )
    ch_alignment = FGUMI_RUNALL.out.bam

    if (!skip_deduplication) {
        FGUMI_DEDUP(ch_alignment)
        ch_alignment = FGUMI_DEDUP.out.bam
        ch_dedup_metrics = FGUMI_DEDUP.out.metrics
    }

    // Coordinate-sort and index
    FGUMI_SORT(ch_alignment)

    BAM_STATS_SAMTOOLS(
        FGUMI_SORT.out.bam.join(FGUMI_SORT.out.index),
        ch_fasta_fai,
    )

    ch_multiqc_files = BAM_STATS_SAMTOOLS.out.stats
        .mix(BAM_STATS_SAMTOOLS.out.flagstat)
        .mix(BAM_STATS_SAMTOOLS.out.idxstats)
        .collect { _meta, file -> file }

    emit:
    bam     = FGUMI_SORT.out.bam   // channel: [ val(meta), bam ]
    index   = FGUMI_SORT.out.index // channel: [ val(meta), bai ]
    metrics = ch_dedup_metrics     // channel: [ val(meta), metrics ]
    multiqc = ch_multiqc_files     // channel: [ files ]
}
