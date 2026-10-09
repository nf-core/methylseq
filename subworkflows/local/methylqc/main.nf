/*
 * EvolvDx methylation QC subworkflow
 *
 * Merges per-sample MethylDackel outputs (methylKit, bedGraph, M-bias) into cohort-level
 * objects, then runs control statistics, read statistics, tiled statistics, Twist and
 * repeat annotation summaries and per-chromosome methylation histograms. The tables are
 * converted into MultiQC custom content (emit: multiqc_files).
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
include { METHYLQC_MULTIQC              } from '../../../modules/local/exo/methylqc_multiqc'

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

    // Collect per-sample files into cohort-level lists, sorted by sample id so that the
    // file list and the metadata list always line up and the task inputs are stable on -resume
    ch_methylkit_cohort = cohort(ch_methylkit)
    ch_bedgraph_cohort  = cohort(ch_bedgraph)
    ch_mbias_cohort     = cohort(ch_mbias)

    COMBINED_METHOBJ (
        ch_methylkit_cohort.files,
        ch_methylkit_cohort.meta
    )
    ch_versions = ch_versions.mix(COMBINED_METHOBJ.out.versions)

    COMBINED_HDF5 (
        ch_bedgraph_cohort.files,
        ch_bedgraph_cohort.meta
    )
    ch_versions = ch_versions.mix(COMBINED_HDF5.out.versions)

    CONTROL_STATS (
        COMBINED_METHOBJ.out.methobj
    )
    ch_versions = ch_versions.mix(CONTROL_STATS.out.versions)

    BASIC_READ_STATISTICS (
        COMBINED_HDF5.out.merged_bedgraph_hdf5,
        ch_mbias_cohort.files,
        ch_mbias_cohort.meta
    )
    ch_versions = ch_versions.mix(BASIC_READ_STATISTICS.out.versions)

    TILE_METHYL_COUNTS (
        COMBINED_METHOBJ.out.methobj_lowcov
    )
    ch_versions = ch_versions.mix(TILE_METHYL_COUNTS.out.versions)

    MERGE_TILED_STATS (
        TILE_METHYL_COUNTS.out.methobj_lowcov_tiled
    )
    ch_versions = ch_versions.mix(MERGE_TILED_STATS.out.versions)

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
        ch_bedgraph_cohort.files,
        ch_bedgraph_cohort.meta,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(METH_HISTOGRAM_BY_CHROMOSOME.out.versions)

    METHYLQC_MULTIQC (
        CONTROL_STATS.out.control_stats,
        BASIC_READ_STATISTICS.out.metrics_table,
        BASIC_READ_STATISTICS.out.meth_distribution,
        BASIC_READ_STATISTICS.out.depth_distribution,
        BASIC_READ_STATISTICS.out.merged_mbias,
        METH_HISTOGRAM_BY_CHROMOSOME.out.chr_percent_df,
        MERGE_ANNOTATION_STATS.out.twist_annotation,
        MERGE_ANNOTATION_STATS.out.repeat_annotation,
        MERGE_TILED_STATS.out.merge_tiled_stats
    )
    ch_versions = ch_versions.mix(METHYLQC_MULTIQC.out.versions)

    emit:
    methobj        = COMBINED_METHOBJ.out.methobj          // channel: path(methobj.rds)
    methobj_lowcov = COMBINED_METHOBJ.out.methobj_lowcov   // channel: path(methobj_lowcov.rds)
    control_stats  = CONTROL_STATS.out.control_stats       // channel: path(control_stats.csv)
    multiqc_files  = METHYLQC_MULTIQC.out.mqc              // channel: path(*_mqc.json)
    versions       = ch_versions                           // channel: path(versions.yml)
}

/*
 * Turn a per-sample channel of [ meta, files ] into two value channels holding the
 * cohort's files and metadata, both ordered by sample id.
 */
def cohort(ch_samples) {
    def ch_sorted = ch_samples.toSortedList { a, b -> a[0].id <=> b[0].id }
    return [
        files: ch_sorted.map { samples -> samples.collect { _meta, files -> files }.flatten() },
        meta : ch_sorted.map { samples -> samples.collect { meta, _files -> meta } }
    ]
}
