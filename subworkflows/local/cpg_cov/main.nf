/*
 * CpG island coverage subworkflow
 *
 * Intersects the alignments with a CpG island BED and reports per-base coverage.
 */

include { BEDTOOLS_INTERSECT_CPG     } from '../../../modules/local/exo/bedtools/intersect_cpg/main'
include { BEDTOOLS_PERBASE_GENOMECOV } from '../../../modules/local/exo/bedtools/per_base_genomecov/main'

workflow CPG_COV {

    take:
    ch_alignments   // channel: [ val(meta), path(bam) ]
    ch_cpg_bed      // channel: path(cpg_island.bed) (value channel)

    main:
    ch_versions = channel.empty()

    BEDTOOLS_INTERSECT_CPG (
        ch_alignments,
        'bam',
        ch_cpg_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_CPG.out.version.first())

    BEDTOOLS_PERBASE_GENOMECOV (
        BEDTOOLS_INTERSECT_CPG.out.intersect,
        'txt'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_PERBASE_GENOMECOV.out.version.first())

    emit:
    intersect = BEDTOOLS_INTERSECT_CPG.out.intersect      // channel: [ val(meta), path(bam) ]
    genomecov = BEDTOOLS_PERBASE_GENOMECOV.out.genomecov  // channel: [ val(meta), path(txt) ]
    versions  = ch_versions                               // channel: path(versions.yml)
}
