/*
 * targeted_sequencing subworkflow

   Filters bedGraph files with the target_regions BED file so they contain positions only listed in the BED.
   If specified, it also generates some performance metrics with Picard CollectHsMetrics, such as Fold-80 Base Penalty,
   HS Library Size, Percent Duplicates, and Percent Off Bait. This is relevant for methylome experiments with targeted seq.
 */

include { BEDTOOLS_INTERSECT              } from '../../../modules/nf-core/bedtools/intersect/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_BEDTOINTERVALLIST        } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS         } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { SAMTOOLS_INDEX                  } from '../../../modules/nf-core/samtools/index/main'
include { ensureList                      } from '../utils_nfcore_methylseq_pipeline'


workflow TARGETED_SEQUENCING {
    
    take:
    ch_bedgraph            // channel: [ val(meta), [ bedGraph(s) ]] when bwameth, [ val(meta), bedGraph ] when bismark
    target_regions        
    ch_fasta               // channel: [ [:], /path/to/genome.fa]
    ch_fasta_index         // channel: /path/to/genome.fa.fai
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment 
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment

    main:

    /*
    * Initial checks
    */
    // for the targeted analysis, the --target_regions_file must be defines
    if (target_regions == null){
        error "Error: --target_regions_file is required when --run_targeted_sequencing is set to true."
    }
    else {
        // check if the --target_regions_file exist
        if (!file(target_regions).exists()) {
            error "Error: The specified region bed file '${target_regions}' does not exist."
        }
    }
    versions = Channel.empty()

    /*
    * Intersect bedGraph files with target regions
    */
    // ensure ch_bedgraph_array contains the bedGraph file(s) in an array
    ch_bedgraph_array = ch_bedgraph
      .map { meta, bedgrahps -> tuple(meta, ensureList(bedgrahps)) }
    // split bedGraphs in separate channel and add the target_regions file
    ch_bedgraphs_target = ch_bedgraph_array
      .flatMap { meta, bedgraphs -> bedgraphs.collect{ bedgraph -> [meta, bedgraph] }}
      .combine(Channel.fromPath(params.target_regions_file))
         
    BEDTOOLS_INTERSECT(ch_bedgraphs_target, [[:], []])
    versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)

    /*
    * Run Picard CollectHSMetrics
    */
    picard_metrics = Channel.empty()
    if (params.run_picard_collecthsmetrics){
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
        picard_metrics = PICARD_COLLECTHSMETRICS.out.metrics

        versions = versions
          .mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
          .mix(PICARD_BEDTOINTERVALLIST.out.versions)
          .mix(PICARD_COLLECTHSMETRICS.out.versions)
    }

    emit:
    picard_metrics        // tuple val(meta), path("*_metrics")
    versions              // path: *.version.txt
}