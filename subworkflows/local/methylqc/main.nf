/*
 * EvolvDx methylation QC subworkflow
 *
 * Merges per-sample MethylDackel outputs (methylKit, bedGraph, M-bias) into cohort-level
 * objects, then runs control statistics, read statistics, tiled statistics, Twist and
 * repeat annotation summaries and per-chromosome methylation histograms.
 */

include { COMBINED_METHOBJ              } from '../../../modules/local/exo/make_combined_methobj'
include { COMBINED_HDF5                 } from '../../../modules/local/exo/make_combined_hdf5'
include { CONTROL_STATS                 } from '../../../modules/local/exo/control_stats'
include { BASIC_READ_STATISTICS         } from '../../../modules/local/exo/basic_read_statistics'
include { TILE_METHYL_COUNTS            } from '../../../modules/local/exo/tile_methyl_counts'
include { MERGE_TILED_STATS             } from '../../../modules/local/exo/merge_tiled_stats'
include { TWIST_ANNOTATION              } from '../../../modules/local/exo/twist_annotation'
include { REPEAT_ANNOTATION             } from '../../../modules/local/exo/repeat_annotation'
include { MERGE_ANNOTATION_STATS        } from '../../../modules/local/exo/merge_annotation_stats'
include { METH_HISTOGRAM_BY_CHROMOSOME  } from '../../../modules/local/exo/meth_histogram_by_chromosome'

workflow METHYLQC {

    take:
    ch_methylkit          // channel: [ val(meta), path(methylkit) ]
    ch_bedgraph           // channel: [ val(meta), path(bedgraph) ]
    ch_mbias              // channel: [ val(meta), path(mbias) ]
    ch_twist_bed          // channel: path(twist_methylome.bed)   (value channel)
    ch_repeat_annot       // channel: path(repeat_annotation)     (value channel)
    ch_chrom_sizes        // channel: path(chrom.sizes)           (value channel)

    main:
    ch_versions = channel.empty()

    // Collect per-sample files and sample metadata into cohort-level lists
    ch_methylkit_files = ch_methylkit.collect { _meta, files -> files }
    ch_methylkit_meta  = ch_methylkit.collect { meta, _files -> meta }
    ch_bedgraph_files  = ch_bedgraph.collect { _meta, files -> files }
    ch_bedgraph_meta   = ch_bedgraph.collect { meta, _files -> meta }
    ch_mbias_files     = ch_mbias.collect { _meta, files -> files }
    ch_mbias_meta      = ch_mbias.collect { meta, _files -> meta }

    COMBINED_METHOBJ (
        ch_methylkit_files,
        ch_methylkit_meta
    )
    ch_versions = ch_versions.mix(COMBINED_METHOBJ.out.versions)

    COMBINED_HDF5 (
        ch_bedgraph_files,
        ch_bedgraph_meta
    )
    ch_versions = ch_versions.mix(COMBINED_HDF5.out.versions)

    CONTROL_STATS (
        COMBINED_METHOBJ.out.methobj
    )
    ch_versions = ch_versions.mix(CONTROL_STATS.out.versions)

    BASIC_READ_STATISTICS (
        COMBINED_HDF5.out.merged_bedgraph_hdf5,
        ch_mbias_files,
        ch_mbias_meta
    )
    ch_versions = ch_versions.mix(BASIC_READ_STATISTICS.out.versions)

    TILE_METHYL_COUNTS (
        COMBINED_METHOBJ.out.methobj_lowcov
    )
    ch_versions = ch_versions.mix(TILE_METHYL_COUNTS.out.versions)

    MERGE_TILED_STATS (
        TILE_METHYL_COUNTS.out.methobj_lowcov_tiled,
        ch_methylkit_meta
    )

    TWIST_ANNOTATION (
        COMBINED_METHOBJ.out.methobj_lowcov,
        ch_twist_bed
    )
    ch_versions = ch_versions.mix(TWIST_ANNOTATION.out.versions)

    REPEAT_ANNOTATION (
        COMBINED_METHOBJ.out.methobj_lowcov,
        ch_repeat_annot
    )
    ch_versions = ch_versions.mix(REPEAT_ANNOTATION.out.versions)

    MERGE_ANNOTATION_STATS (
        COMBINED_METHOBJ.out.methobj_lowcov,
        TWIST_ANNOTATION.out.methobj_lowcov_twist_annotation_obj,
        TWIST_ANNOTATION.out.cpg_obj_twist_out,
        REPEAT_ANNOTATION.out.methobj_lowcov_repeat_annotation_obj
    )
    ch_versions = ch_versions.mix(MERGE_ANNOTATION_STATS.out.versions)

    METH_HISTOGRAM_BY_CHROMOSOME (
        ch_bedgraph_files,
        ch_bedgraph_meta,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(METH_HISTOGRAM_BY_CHROMOSOME.out.versions)

    emit:
    methobj        = COMBINED_METHOBJ.out.methobj          // channel: path(methobj.rds)
    methobj_lowcov = COMBINED_METHOBJ.out.methobj_lowcov   // channel: path(methobj_lowcov.rds)
    control_stats  = CONTROL_STATS.out.control_stats       // channel: path(control_stats.csv)
    versions       = ch_versions                           // channel: path(versions.yml)
}
