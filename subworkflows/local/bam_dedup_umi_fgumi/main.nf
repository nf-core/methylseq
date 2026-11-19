//
// UMI-aware duplicate marking with fgumi: copy the UMIs from the unmapped BAM onto the alignments (zipper),
// sort by template coordinate, mark duplicates per UMI family and strand, then coordinate-sort and index
//

include { SAMTOOLS_DICT                         } from '../../../modules/nf-core/samtools/dict/main'
include { FGUMI_SORT as FGUMI_SORT_ALIGNED      } from '../../../modules/nf-core/fgumi/sort/main'
include { FGUMI_SORT as FGUMI_SORT_UNMAPPED     } from '../../../modules/nf-core/fgumi/sort/main'
include { FGUMI_ZIPPER                          } from '../../../modules/nf-core/fgumi/zipper/main'
include { FGUMI_SORT as FGUMI_SORT_TEMPLATE     } from '../../../modules/nf-core/fgumi/sort/main'
include { FGUMI_DEDUP                           } from '../../../modules/nf-core/fgumi/dedup/main'
include { FGUMI_SORT as FGUMI_SORT_COORDINATE   } from '../../../modules/nf-core/fgumi/sort/main'

workflow BAM_DEDUP_UMI_FGUMI {
    take:
    ch_bam       // channel: [ val(meta), bam ], aligned
    ch_ubam      // channel: [ val(meta), bam ], unmapped, from fgumi extract
    ch_fasta_fai // channel: [ val(meta), fasta, fai ]

    main:
    SAMTOOLS_DICT(ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] })

    // zipper needs both BAMs in the same read name order; the aligned BAM arrives coordinate-sorted
    FGUMI_SORT_ALIGNED(ch_bam)
    FGUMI_SORT_UNMAPPED(ch_ubam)

    FGUMI_ZIPPER(
        FGUMI_SORT_ALIGNED.out.bam.join(FGUMI_SORT_UNMAPPED.out.bam, failOnDuplicate: true, failOnMismatch: true),
        ch_fasta_fai.join(SAMTOOLS_DICT.out.dict).first(),
    )

    FGUMI_SORT_TEMPLATE(FGUMI_ZIPPER.out.bam)

    FGUMI_DEDUP(FGUMI_SORT_TEMPLATE.out.bam)

    FGUMI_SORT_COORDINATE(FGUMI_DEDUP.out.bam)

    emit:
    bam       = FGUMI_SORT_COORDINATE.out.bam   // channel: [ val(meta), bam ]
    index     = FGUMI_SORT_COORDINATE.out.index // channel: [ val(meta), bai ]
    metrics   = FGUMI_DEDUP.out.metrics         // channel: [ val(meta), metrics ]
    histogram = FGUMI_DEDUP.out.histogram       // channel: [ val(meta), histogram ]
}
