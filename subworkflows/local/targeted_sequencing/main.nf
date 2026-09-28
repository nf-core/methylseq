/*
 * targeted_sequencing subworkflow
 *
 * Filters bedGraph and coverage files with the target_regions BED file so they contain positions only listed in the BED.
 * If specified, it also generates some performance metrics with Picard CollectHsMetrics, such as Fold-80 Base Penalty,
 * HS Library Size, Percent Duplicates, and Percent Off Bait. This is relevant for methylome experiments with targeted seq.
 */

nextflow.enable.types = true

include { FILTER_BEDGRAPH_TARGETS                      } from '../../../modules/local/filter_bedgraph_targets/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_COV } from '../../../modules/nf-core/bedtools/intersect/main'
include { PICARD_CREATESEQUENCEDICTIONARY              } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_BEDTOINTERVALLIST                     } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS                      } from '../../../modules/nf-core/picard/collecthsmetrics/main'

workflow TARGETED_SEQUENCING {
    take:
    ch_inputs: Channel<TargetedSequencingInput>
    val_target_regions: Value<Path>
    val_fasta: Value<Path>
    val_fasta_index: Value<Path>
    collecthsmetrics: Boolean // whether to run Picard CollectHsMetrics

    main:

    /*
     * Intersect bedGraph files with target regions (CpG-aware boundary handling)
     * Split the bedGraph file(s) of each sample into individual bedGraphs.
     * The FILTER_BEDGRAPH_TARGETS process extends single-C intervals by 1 bp before
     * intersection so that CpGs straddling a target boundary are not lost.
     */
    ch_bedgraphs_target = ch_inputs
        .flatMap { r -> r.bedgraphs.collect { bedgraph -> record(id: r.id, meta: r.meta, bedgraph: bedgraph) } }
        .combine(targets: val_target_regions)

    ch_bedgraph_intersect = FILTER_BEDGRAPH_TARGETS(ch_bedgraphs_target)
        .map { r -> tuple(r.id, r.bedgraph_intersect) }
        .groupBy()
        .map { id, bedgraphs -> record(id: id, bedgraph_intersect: bedgraphs) }

    /*
     * Intersect Bismark coverage files with target regions
     * The .cov.gz files are filtered the same way as bedGraph files so that
     * downstream tools (methylKit, bsseq, DSS) receive only on-target CpGs.
     */
    ch_coverage_target = ch_inputs
        .filter { r -> r.coverage != null }
        .map { r -> record(id: r.id, meta: r.meta, intervals1: r.coverage) }
        .combine(intervals2: val_target_regions)

    ch_coverage_intersect = BEDTOOLS_INTERSECT_COV(ch_coverage_target, null)

    ch_results = ch_bedgraph_intersect
        .join(ch_coverage_intersect, by: 'id', remainder: true)

    /*
     * Run Picard CollectHSMetrics
     */
    if (collecthsmetrics) {
        /*
         * Creation of a dictionary for the reference genome
         */
        val_reference_dict = PICARD_CREATESEQUENCEDICTIONARY(val_fasta)

        /*
         * Conversion of the covered targets BED file to an interval list
         */
        val_intervallist = PICARD_BEDTOINTERVALLIST(val_target_regions, val_reference_dict)

        /*
         * Generation of the metrics
         * Note: Using the same intervals for both target and bait as they are typically
         * the same for targeted methylation sequencing experiments
         */
        ch_picard_hsmetrics = PICARD_COLLECTHSMETRICS(
            ch_inputs.combine(
                bait_intervals: val_intervallist,
                target_intervals: val_intervallist,
                ref: val_fasta,
                ref_fai: val_fasta_index,
                ref_dict: val_reference_dict
            )
        )
        ch_results = ch_results.join(ch_picard_hsmetrics, by: 'id')
    }
    else {
        val_reference_dict = null
        val_intervallist = null
    }

    emit:
    results        : Channel<TargetedSequencingResult> = ch_results
    reference_dict : Value<Path>? = val_reference_dict
    intervallist   : Value<Path>? = val_intervallist
}

record TargetedSequencingInput {
    id: String
    meta: Record
    bam: Path
    bai: Path
    bedgraphs: List<Path>
    coverage: Path?
}

record TargetedSequencingResult {
    id: String
    bedgraph_intersect: Bag<Path>
    coverage_intersect: Path?
    picard_hsmetrics: Path?
}
