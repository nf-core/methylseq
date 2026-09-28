nextflow.enable.types = true

include { METHYLDACKEL_EXTRACT } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS   } from '../../../modules/nf-core/methyldackel/mbias/main'

include { AlignedSample } from '../../../utils/types.nf'

workflow BAM_METHYLDACKEL {

    take:
    ch_bam: Channel<AlignedSample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>

    main:

    /*
     * Extract per-base methylation and plot methylation bias
     */
    ch_inputs = ch_bam.combine(fasta: val_fasta, fai: val_fasta_index)

    ch_extract = METHYLDACKEL_EXTRACT(ch_inputs)

    ch_mbias = METHYLDACKEL_MBIAS(ch_inputs)

    ch_results = ch_extract.join(ch_mbias, by: 'id')

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_results.flatMap { r -> r.methyldackel_bedgraph + r.methyldackel_methylkit + [r.methyldackel_mbias].toSet() }

    emit:
    results : Channel<MethyldackelResult> = ch_results
    multiqc : Channel<Path>               = ch_multiqc_files
}

record MethyldackelResult {
    id: String
    methyldackel_bedgraph: Set<Path>
    methyldackel_methylkit: Set<Path>
    methyldackel_mbias: Path
}
