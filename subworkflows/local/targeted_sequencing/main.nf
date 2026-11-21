/*
 * targeted_sequencing subworkflow
 *
 * Filters bedGraph and coverage files to the target_regions BED file, keeping both strands of a CpG cut by a target edge.
 * If specified, it also generates some performance metrics with Picard CollectHsMetrics, such as Fold-80 Base Penalty,
 * HS Library Size, Percent Duplicates, and Percent Off Bait. This is relevant for methylome experiments with targeted seq.
 */

include { WIDEN_TARGETS_CPG                            } from '../../../modules/local/widen_targets_cpg/main'
include { FILTER_BEDGRAPH_TARGETS                      } from '../../../modules/local/filter_bedgraph_targets/main'
include { FILTER_COV_TARGETS                           } from '../../../modules/local/filter_cov_targets/main'
include { PICARD_CREATESEQUENCEDICTIONARY              } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_BEDTOINTERVALLIST                     } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS                      } from '../../../modules/nf-core/picard/collecthsmetrics/main'

workflow TARGETED_SEQUENCING {
    take:
    ch_bedgraph // channel: [ val(meta), [ bedGraph(s) ]] when bwameth, [ val(meta), bedGraph ] when bismark
    ch_coverage // channel: [ val(meta), cov.gz ] from Bismark methylation extractor (empty for bwameth)
    ch_target_regions // channel: path(target_regions.bed)
    ch_fasta // channel: [ [:], /path/to/genome.fa]
    ch_fasta_index // channel: [ val(meta), /path/to/genome.fa.fai]
    ch_bam // channel: [ val(meta), [ bam ] ] ## BAM from alignment
    ch_bai // channel: [ val(meta), [ bai ] ] ## BAI from alignment
    ch_gzi // channel: [ val(meta), [ gzi ] ] ## GZI from fasta
    collecthsmetrics // boolean: whether to run Picard CollectHsMetrics
    taps // boolean: TAPS run (ch_bedgraph holds rastair CpG calls)

    main:

    ch_versions = Channel.empty()
    ch_picard_metrics = Channel.empty()

    // rastair calls carry their strand, so TAPS runs filter against the targets as given
    ch_filter_targets = ch_target_regions
    if (!taps) {
        WIDEN_TARGETS_CPG(
            ch_target_regions,
            ch_fasta.ifEmpty { error("The targeted analysis needs the FASTA to widen the targets at cut CpGs") },
            ch_fasta_index.ifEmpty { error("The targeted analysis needs the FASTA index to widen the targets at cut CpGs") },
        )
        ch_filter_targets = WIDEN_TARGETS_CPG.out.targets
        ch_versions = ch_versions.mix(WIDEN_TARGETS_CPG.out.versions)

        WIDEN_TARGETS_CPG.out.issues
            .filter { issues -> issues.size() > 0 }
            .subscribe { issues ->
                def lines = issues.readLines()
                def contigs = lines.collect { line -> line.split('\t', -1)[1] }.unique()
                log.warn(
                    "--target_regions_file: ${lines.size()} target(s) skipped or clipped to fit the FASTA, on contig(s) ${contigs.take(10).join(', ')}${contigs.size() > 10 ? ', ...' : ''}. First ${Math.min(5, lines.size())}:\n" + lines.take(5).join('\n')
                )
            }
    }

    // Split bwameth's list of bedGraphs into one task per file
    ch_bedgraphs_target = ch_bedgraph
        .map { meta, bedgraphs -> tuple(meta, bedgraphs instanceof List ? bedgraphs : [bedgraphs]) }
        .flatMap { meta, bedgraphs -> bedgraphs.collect { bedgraph -> [meta, bedgraph] } }
        .combine(ch_filter_targets)

    FILTER_BEDGRAPH_TARGETS(ch_bedgraphs_target)
    ch_versions = ch_versions.mix(FILTER_BEDGRAPH_TARGETS.out.versions)

    // Bismark coverage files, for methylKit, bsseq and DSS
    FILTER_COV_TARGETS(ch_coverage.combine(ch_filter_targets))
    ch_versions = ch_versions.mix(FILTER_COV_TARGETS.out.versions)

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
            [],
        )
        ch_intervals = PICARD_BEDTOINTERVALLIST.out.intervallist.map { it[1] }
        ch_versions = ch_versions.mix(PICARD_BEDTOINTERVALLIST.out.versions)

        /*
         * Generation of the metrics
         * Note: Using the same intervals for both target and bait as they are typically
         * the same for targeted methylation sequencing experiments
         */
        ch_picard_inputs = ch_bam
            .join(ch_bai)
            .combine(ch_intervals)
            .combine(ch_intervals)
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_sequence_dictionary)
            .combine(ch_gzi)
            .multiMap { meta, bam, bai, intervals1, intervals2, meta_fasta, fasta, meta_fasta_index, fasta_index, meta_dict, dict, meta_gzi, gzi ->
                bam_etc: [meta, bam, bai, intervals1, intervals2]
                fasta: [meta_fasta, fasta]
                fasta_index: [meta_fasta_index, fasta_index]
                dict: [meta_dict, dict]
                gzi: [meta_gzi, gzi]
            }

        PICARD_COLLECTHSMETRICS(
            ch_picard_inputs.bam_etc,
            ch_picard_inputs.fasta,
            ch_picard_inputs.fasta_index,
            ch_picard_inputs.dict,
            ch_picard_inputs.gzi,
        )
        ch_picard_metrics = PICARD_COLLECTHSMETRICS.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions)
    }

    emit:
    bedgraph_filtered = FILTER_BEDGRAPH_TARGETS.out.intersect // channel: [ val(meta), path("*.targeted.{bedGraph,bed}") ]
    coverage_filtered = FILTER_COV_TARGETS.out.intersect // channel: [ val(meta), path("*.targeted.cov") ]
    picard_metrics    = ch_picard_metrics // channel: [ val(meta), path("*_metrics") ]
    versions          = ch_versions // channel: path("*.version.txt")
}
