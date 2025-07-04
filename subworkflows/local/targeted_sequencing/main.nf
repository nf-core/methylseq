/*
 * targeted_sequencing subworkflow
 *
 * Filters bedGraph files with the target_regions BED file so they contain positions only listed in the BED.
 * If specified, it also generates some performance metrics with Picard CollectHsMetrics, such as Fold-80 Base Penalty,
 * HS Library Size, Percent Duplicates, and Percent Off Bait. This is relevant for methylome experiments with targeted seq.
 */

include { BEDTOOLS_INTERSECT              } from '../../../modules/nf-core/bedtools/intersect/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_BEDTOINTERVALLIST        } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS         } from '../../../modules/nf-core/picard/collecthsmetrics/main'

workflow TARGETED_SEQUENCING {

    take:
    ch_bedgraph            // channel: [ val(meta), [ bedGraph(s) ]] when bwameth, [ val(meta), bedGraph ] when bismark
    ch_target_regions      // channel: path(target_regions.bed)
    ch_fasta               // channel: [ [:], /path/to/genome.fa]
    ch_fasta_index         // channel: [ val(meta), /path/to/genome.fa.fai]
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment
    collecthsmetrics       // boolean: whether to run Picard CollectHsMetrics

    main:

    ch_versions = Channel.empty()
    ch_picard_metrics = Channel.empty()

    /*
     * Intersect bedGraph files with target regions
     * Ensure ch_bedgraph contains the bedGraph file(s) in an array and split into individual bedGraphs
     */
    ch_bedgraphs_target = ch_bedgraph
        .map { meta, bedgraphs -> tuple(meta, bedgraphs instanceof List ? bedgraphs : [bedgraphs]) }
        .flatMap { meta, bedgraphs -> bedgraphs.collect { bedgraph -> [meta, bedgraph] } }
        .combine(ch_target_regions)

    BEDTOOLS_INTERSECT(
        ch_bedgraphs_target,
        [[:], []]
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

    /*
     * Run Picard CollectHSMetrics
     */
    if (collecthsmetrics) {
        // Create target regions with meta for Picard tools
        target_regions_with_meta = ch_target_regions.map { target_file ->
            tuple(["id": file(target_file).baseName], target_file)
        }

        /*
         * Creation of a dictionary for the reference genome
         */
        PICARD_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_sequence_dictionary = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
        ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

        /*
         * Conversion of the covered targets BED file to an interval list
         */
        PICARD_BEDTOINTERVALLIST(
            target_regions_with_meta,
            ch_sequence_dictionary,
            []
        )
        ch_intervals = PICARD_BEDTOINTERVALLIST.out.intervallist.map { it[1] }
        ch_versions = ch_versions.mix(PICARD_BEDTOINTERVALLIST.out.versions)

        /*
         * Generation of the metrics
         * Note: Using the same intervals for both target and bait as they are typically
         * the same for targeted methylation sequencing experiments
         */
        ch_picard_inputs = ch_bam.join(ch_bai)
            .combine(ch_intervals)
            .combine(ch_intervals)
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_sequence_dictionary)
            .multiMap { meta, bam, bai, intervals1, intervals2, meta_fasta, fasta, meta_fasta_index, fasta_index, meta_dict, dict ->
                bam_etc: [ meta, bam, bai, intervals1, intervals2 ] // intervals: baits, targets
                fasta: [ meta_fasta, fasta ]
                fasta_index: [ meta_fasta_index, fasta_index ]
                dict: [ meta_dict, dict ]
            }

        PICARD_COLLECTHSMETRICS(
            ch_picard_inputs.bam_etc,
            ch_picard_inputs.fasta,
            ch_picard_inputs.fasta_index,
            ch_picard_inputs.dict
        )
        ch_picard_metrics = PICARD_COLLECTHSMETRICS.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions)
    }

    emit:
    bedgraph_filtered = BEDTOOLS_INTERSECT.out.intersect  // channel: [ val(meta), path("*.bedGraph") ]
    picard_metrics    = ch_picard_metrics                 // channel: [ val(meta), path("*_metrics") ]
    versions          = ch_versions                       // channel: path("*.version.txt")
}
