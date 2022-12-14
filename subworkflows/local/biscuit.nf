/*
 * biscuit subworkflow
 */

include { BISCUIT_INDEX   } from '../../modules/nf-core/biscuit/index/main'
include { BISCUIT_ALIGN   } from '../../modules/nf-core/biscuit/align/main'
include { BISCUIT_BLASTER } from '../../modules/nf-core/biscuit/biscuitblaster/main'
include { BISCUIT_PILEUP  } from '../../modules/nf-core/biscuit/pileup/main'
include { BISCUIT_QC      } from '../../modules/nf-core/biscuit/qc/main'

workflow BISCUIT {
    take:
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    versions = Channel.empty()

    /*
     * Generate biscuit index if not supplied
     */
    if (params.biscuit_index) {
        biscuit_index = file(params.biscuit_index)
    } else {
        BISCUIT_INDEX(params.fasta)
        biscuit_index = BISCUIT_INDEX.out.index
        versions = versions.mix(BISCUIT_INDEX.out.versions)
    }

    /*
     * Align with biscuit; mark duplicates unless params.skip_deduplication
     */
    if (params.skip_deduplication || params.rrbs){
        BISCUIT_ALIGN (
            reads,
            biscuit_index
        )
        versions = versions.mix(BISCUIT_ALIGN.out.versions)
        alignments = BISCUIT_ALIGN.out.indexed_bam
    } else {
        BISCUIT_BLASTER (
            reads,
            biscuit_index
        )
        versions = versions.mix(BISCUIT_BLASTER.out.versions)
        alignments = BISCUIT_BLASTER.out.indexed_bam
    }

    /*
     * Extract methylation calls
     */
    BISCUIT_PILEUP (
        alignments.map {
            meta, bam, bai -> [ meta, bam, bai, [], [] ] // Supply empty lists for tumor bam/bai (optional inputs to pileup process)
        },
        biscuit_index
    )
    versions = versions.mix(BISCUIT_PILEUP.out.versions)

    /*
     * QC alignments
     */
    // BISCUIT_QC (
        // alignments.map {
        //     meta, bam, bai -> [ meta, bam ] // remove bai index
        // },
        // biscuit_index
    // )

    multiqc_files = Channel.empty()

    emit:
    bam   = alignments
    dedup = alignments
    mqc   = multiqc_files                        // path: *{html,txt}
    versions
}


// TODO:
// methylation extraction
// samtools flagstat, stats
// reports
// multiqc




