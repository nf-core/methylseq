/*
 * targeted_analysis subworkflow

   Generation of some performance metrics with Picard HsMetrics, such as Fold-80 Base Penalty,
   HS Library Size, Percent Duplicates, and Percent Off Bait. This is relevant for methylome 
   experiments with targeted seq.
 */

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_BEDTOINTERVALLIST        } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS         } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { SAMTOOLS_INDEX                  } from '../../../modules/nf-core/samtools/index/main'

workflow PICARD_TARGETED_SEQUENCING {
    
    take:
    ch_fasta               // channel: /path/to/genome.fa
    ch_fasta_index         // channel: /path/to/genome.fa.fai
    target_regions         // channel: /path/to/target_regions.bed ## BED file with the covered targets
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment 
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment

    main:
    versions = Channel.empty()

    targets_with_meta = tuple(["id": file(target_regions).name.replaceFirst(~/\.[^\.]+$/, '')], target_regions)

    /*
     * Creation of a dictionary for the reference genome
     */
    PICARD_CREATESEQUENCEDICTIONARY(ch_fasta)
    versions = versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    /*
     * Conversion of the covered targets BED file to an interval list
     */
    PICARD_BEDTOINTERVALLIST (
        targets_with_meta,
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
        []
    )
    versions = versions.mix(PICARD_BEDTOINTERVALLIST.out.versions)
    interval_file = PICARD_BEDTOINTERVALLIST.out.interval_list.map{ it[1] }
    
    /*
     * Generation of the metrics
     */
    PICARD_COLLECTHSMETRICS (
       ch_bam.join(ch_bai),
       ch_fasta,
       ch_fasta_index,
       interval_file
    )
    versions = versions.mix(PICARD_COLLECTHSMETRICS.out.versions)

    emit:
    metrics  = PICARD_COLLECTHSMETRICS.out.metrics    // tuple val(meta), path("*_metrics")
    versions                                          // path: *.version.txt

}