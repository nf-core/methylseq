//
// Move UMIs out of the reads with fgumi: FASTQ -> unmapped BAM (UMI in the RX tag) -> FASTQ without the UMI bases
//

include { FGUMI_EXTRACT } from '../../../modules/nf-core/fgumi/extract/main'
include { FGUMI_FASTQ   } from '../../../modules/nf-core/fgumi/fastq/main'

workflow FASTQ_EXTRACT_UMI_FGUMI {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    FGUMI_EXTRACT(
        ch_reads.map { meta, reads -> [meta, reads, meta.id] }
    )

    FGUMI_FASTQ(FGUMI_EXTRACT.out.bam)

    // Paired reads are written with --out1/--out2 (see conf/modules/fgumi_fastq.config);
    // the module still gzips fgumi's then-empty stdout stream, so drop that file
    ch_umi_reads = FGUMI_FASTQ.out.reads.map { meta, reads ->
        def files = reads instanceof List ? reads : [reads]
        [meta, meta.single_end ? files : files.findAll { fastq -> fastq.name ==~ /.*_[12]\.fastq\.gz/ }.sort { fastq -> fastq.name }]
    }

    emit:
    reads = ch_umi_reads           // channel: [ val(meta), [ reads ] ]
    ubam  = FGUMI_EXTRACT.out.bam  // channel: [ val(meta), unmapped bam ]
}
