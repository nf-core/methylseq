//
// UMI-aware duplicate marking with fgumi: copy the UMIs from the unmapped BAM onto the alignments (zipper),
// sort by template coordinate, mark duplicates per UMI family and strand, then coordinate-sort, index and run samtools stats
//

include { SAMTOOLS_DICT                         } from '../../../modules/nf-core/samtools/dict/main'
include { FGUMI_ZIPPER                          } from '../../../modules/nf-core/fgumi/zipper/main'
include { FGUMI_SORT as FGUMI_SORT_TEMPLATE     } from '../../../modules/nf-core/fgumi/sort/main'
include { FGUMI_DEDUP                           } from '../../../modules/nf-core/fgumi/dedup/main'
include { FGUMI_SORT as FGUMI_SORT_COORDINATE   } from '../../../modules/nf-core/fgumi/sort/main'
include { BAM_STATS_SAMTOOLS                    } from '../../nf-core/bam_stats_samtools/main'

workflow BAM_DEDUP_UMI_FGUMI {
    take:
    ch_bam       // channel: [ val(meta), bam ], aligned, unsorted (input read order)
    ch_ubam      // channel: [ val(meta), bam ], unmapped, from fgumi extract
    ch_fasta_fai // channel: [ val(meta), fasta, fai ]

    main:
    SAMTOOLS_DICT(ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] })

    // zipper needs both BAMs in the same read name order. Both are in input read order, as Trim Galore keeps it,
    // so neither is sorted first
    FGUMI_ZIPPER(
        ch_bam.join(ch_ubam, failOnDuplicate: true, failOnMismatch: true),
        ch_fasta_fai.join(SAMTOOLS_DICT.out.dict).first(),
    )

    FGUMI_SORT_TEMPLATE(FGUMI_ZIPPER.out.bam)

    FGUMI_DEDUP(FGUMI_SORT_TEMPLATE.out.bam)

    FGUMI_SORT_COORDINATE(FGUMI_DEDUP.out.bam)

    BAM_STATS_SAMTOOLS(
        FGUMI_SORT_COORDINATE.out.bam.join(FGUMI_SORT_COORDINATE.out.index),
        ch_fasta_fai.first(),
    )

    ch_multiqc_files = BAM_STATS_SAMTOOLS.out.stats
        .mix(BAM_STATS_SAMTOOLS.out.flagstat)
        .mix(BAM_STATS_SAMTOOLS.out.idxstats)
        .collect { _meta, file -> file }

    emit:
    bam       = FGUMI_SORT_COORDINATE.out.bam   // channel: [ val(meta), bam ]
    index     = FGUMI_SORT_COORDINATE.out.index // channel: [ val(meta), bai ]
    metrics   = FGUMI_DEDUP.out.metrics         // channel: [ val(meta), metrics ]
    histogram = FGUMI_DEDUP.out.histogram       // channel: [ val(meta), histogram ]
    multiqc   = ch_multiqc_files                // channel: [ files ]
}
