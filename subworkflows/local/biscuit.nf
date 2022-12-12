/*
 * biscuit subworkflow
 */

include { BISCUIT_INDEX } from '../../modules/nf-core/biscuit/index/main'
include { BISCUIT_BLASTER } from '../../modules/nf-core/biscuit/biscuitblaster/main'

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
     * Align with biscuit
     */
    BISCUIT_BLASTER (
        reads,
        biscuit_index
    )
    versions = versions.mix(BISCUIT_BLASTER.out.versions)

    multiqc_files = Channel.empty()

    emit:
    bam   = BISCUIT_BLASTER.out.bam
    dedup = BISCUIT_BLASTER.out.bam
    mqc   = multiqc_files                        // path: *{html,txt}
    versions
}
