/*
 * targeted_sequencing subworkflow
 *
 * Filters bedGraph files with the target_regions BED file so they contain positions only listed in the BED.
 * If specified, it also generates some performance metrics with Picard CollectHsMetrics, such as Fold-80 Base Penalty,
 * HS Library Size, Percent Duplicates, and Percent Off Bait. This is relevant for methylome experiments with targeted seq.
 */

nextflow.preview.types = true

include { BEDTOOLS_INTERSECT              } from '../../../modules/nf-core/bedtools/intersect/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_BEDTOINTERVALLIST        } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS         } from '../../../modules/nf-core/picard/collecthsmetrics/main'

workflow TARGETED_SEQUENCING {

    take:
    ch_inputs: Channel<AlignedSample>
    val_target_regions: Value<Path>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>
    collecthsmetrics: Boolean

    main:

    /*
     * Intersect bedGraph files with target regions
     * Ensure ch_bedgraph contains the bedGraph file(s) in an array and split into individual bedGraphs
     */
    ch_intersect_inputs = ch_inputs
        .map { r ->
            record(
                id: r.id,
                intervals1: r.bedgraph,
            )
        }
        .combine(intervals2: val_target_regions)

    ch_results = BEDTOOLS_INTERSECT( ch_intersect_inputs )

    /*
     * Run Picard CollectHSMetrics
     */
    if (collecthsmetrics) {

        /*
         * Creation of a dictionary for the reference genome
         */
        val_reference_dict_inputs = val_fasta.map { fa -> record(id: fa.baseName, fasta: fa) }
        val_reference_dict = PICARD_CREATESEQUENCEDICTIONARY( val_reference_dict_inputs ).map { r -> r.reference_dict }

        /*
         * Conversion of the covered targets BED file to an interval list
         */
        val_intervallist_inputs = val_target_regions.map { tr -> record(id: tr.baseName, bed: tr) }
        val_intervallist = PICARD_BEDTOINTERVALLIST( val_intervallist_inputs.combine(reference_dict: val_reference_dict) ).map { r -> r.intervallist }

        /*
         * Generation of the metrics
         * Note: Using the same intervals for both target and bait as they are typically
         * the same for targeted methylation sequencing experiments
         */
        ch_picard_hsmetrics = PICARD_COLLECTHSMETRICS(
            ch_inputs.combine(
                bait_intervals: val_intervallist,
                target_intervals: val_intervallist,
                fasta: val_fasta,
                fai: val_fasta_index,
                reference_dict: val_reference_dict
            )
        )
        ch_results = ch_results.join(ch_picard_hsmetrics, by: 'id')
    } else {
        val_reference_dict = null
        val_intervallist = null
    }

    emit:
    results        : Channel<TargetedSequencingResult> = ch_results
    reference_dict : Value<Record>? = val_reference_dict
    intervallist   : Value<Record>? = val_intervallist   
}


record AlignedSample {
    id: String
    bedgraph: Path
    bam: Path
    bai: Path
}

record TargetedSequencingResult {
    id: String
    bedgraph_intersect: Path
    picard_hsmetrics: Path?
}
