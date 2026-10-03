nextflow.enable.types = true

include { METHYLDACKEL_EXTRACT } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS   } from '../../../modules/nf-core/methyldackel/mbias/main'

workflow BAM_METHYLDACKEL {

    take:
    ch_bam: Channel<AlignedSample>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>
    args: MethyldackelArgs

    main:

    /*
     * Extract per-base methylation and plot methylation bias
     */
    ch_inputs = ch_bam.combine(fasta: val_fasta, fai: val_fasta_index)

    ch_extract = METHYLDACKEL_EXTRACT(ch_inputs.map { r -> r + record(args: args.extract) })

    ch_mbias = METHYLDACKEL_MBIAS(ch_inputs.map { r -> r + record(args: args.mbias) })

    ch_results = ch_extract.join(ch_mbias, by: 'meta')

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_results.flatMap { r -> r.methyldackel_bedgraph + r.methyldackel_methylkit + [r.methyldackel_mbias].toSet() }

    emit:
    results : Channel<MethyldackelResult> = ch_results
    multiqc : Channel<Path>               = ch_multiqc_files
}

record MethyldackelArgs {
    extract: String?
    mbias: String?
}

record MethyldackelResult {
    meta: Record
    methyldackel_bedgraph: Set<Path>
    methyldackel_methylkit: Set<Path>
    methyldackel_mbias: Path
}

record AlignedSample {
    meta: Record
    bam: Path
    bai: Path
}
