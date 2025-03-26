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
    ch_fasta               // channel: [ [:], /path/to/genome.fa]
    ch_fasta_index         // channel: /path/to/genome.fa.fai
    target_regions         // channel: /path/to/target_regions.bed ## BED file with the covered targets
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment 
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment

    main:

    // setup channels for versions, fasta, fasta_index and target regions
    versions = Channel.empty()
    ch_fasta_with_meta = ch_fasta
      .map{meta, fasta -> tuple(meta + ["id": file(fasta).name.replaceFirst(~/\.[^\.]+$/, '')], fasta)}
    ch_fasta_index_with_meta = ch_fasta_index.map{fasta_index -> 
      tuple(
        ["id": file(fasta_index).name.replaceFirst(~/\.[^\.]+$/, '')],
        fasta_index
        )
      }
    //fasta_index_with_meta = tuple(["id": file(fasta_index).BaseName], fasta_index)
    target_regions_with_meta = tuple(["id": file(target_regions).BaseName], target_regions)

    /*
     * Creation of a dictionary for the reference genome
     */
    PICARD_CREATESEQUENCEDICTIONARY(ch_fasta_with_meta)
    ch_sequence_dictionary = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict

    /*
     * Conversion of the covered targets BED file to an interval list
     */
    PICARD_BEDTOINTERVALLIST (
        target_regions_with_meta,
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
        []
    )
    ch_intervals = PICARD_BEDTOINTERVALLIST.out.interval_list.map{ it[1] }
    
    /*
     * Generation of the metrics
     */
    PICARD_COLLECTHSMETRICS(
       ch_bam.join(ch_bai).combine(ch_intervals).combine(ch_intervals),
       ch_fasta_with_meta,
       ch_fasta_index_with_meta,
       ch_sequence_dictionary
    )
    versions = versions
      .mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
      .mix(PICARD_BEDTOINTERVALLIST.out.versions)
      .mix(PICARD_COLLECTHSMETRICS.out.versions)

    emit:
    metrics  = PICARD_COLLECTHSMETRICS.out.metrics    // tuple val(meta), path("*_metrics")
    versions                                          // path: *.version.txt

}