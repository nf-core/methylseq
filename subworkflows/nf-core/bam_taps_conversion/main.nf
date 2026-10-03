/*
 * TAPS methylation conversion subworkflow
 *
 * Uses Rastair to assess C->T conversion as a readout for methylation in a genome-wide basis
 */

nextflow.enable.types = true

include { RASTAIR_MBIAS             } from '../../../modules/nf-core/rastair/mbias/main'
include { RASTAIR_MBIASPARSER       } from '../../../modules/nf-core/rastair/mbiasparser/main'
include { RASTAIR_CALL              } from '../../../modules/nf-core/rastair/call/main'
include { RASTAIR_METHYLKIT         } from '../../../modules/nf-core/rastair/methylkit/main'

workflow BAM_TAPS_CONVERSION {

    take:
    ch_bam: Channel<AlignedSample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>

    main:

    log.info "Running TAPS conversion module with Rastair to assess C->T conversion as a readout for methylation."

    ch_inputs = ch_bam.combine(fasta: val_fasta, fai: val_fasta_index)

    ch_rastair_mbias = RASTAIR_MBIAS(ch_inputs)

    ch_rastair_mbiasparser = RASTAIR_MBIASPARSER(ch_rastair_mbias)

    ch_rastair_call = RASTAIR_CALL(ch_inputs.join(ch_rastair_mbiasparser, by: 'meta'))

    ch_rastair_methylkit = RASTAIR_METHYLKIT(ch_rastair_call)

    ch_results = ch_rastair_mbias
        .join(ch_rastair_mbiasparser, by: 'meta')
        .join(ch_rastair_call, by: 'meta')
        .join(ch_rastair_methylkit, by: 'meta')

    emit:
    ch_results
}

record AlignedSample {
    meta: Record
    bam: Path
    bai: Path
}
