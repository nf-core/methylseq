/*
 * biscuit subworkflow
 */

include { BISCUIT_INDEX } from '../../modules/nf-core/biscuit/index/main'
include { BISCUIT_ALIGN } from '../../modules/nf-core/biscuit/align/main'

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
    } else {
        BISCUIT_BLASTER (
            reads,
            biscuit_index
        )
    }

    multiqc_files = Channel.empty()

    emit:
    bam   = BISCUIT_ALIGN.out.bam
    dedup = BISCUIT_ALIGN.out.bam
    mqc   = multiqc_files                        // path: *{html,txt}
    versions
}
