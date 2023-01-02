/*
 * biscuit subworkflow
 */

include { BISCUIT_INDEX                           } from '../../modules/nf-core/biscuit/index/main'
include { BISCUIT_ALIGN                           } from '../../modules/nf-core/biscuit/align/main'
include { BISCUIT_BLASTER                         } from '../../modules/nf-core/biscuit/biscuitblaster/main'
include { BISCUIT_BSCONV                          } from '../../modules/nf-core/biscuit/bsconv/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BSCONV } from '../../modules/nf-core/samtools/index/main'
include { BISCUIT_PILEUP                          } from '../../modules/nf-core/biscuit/pileup/main'
include { BISCUIT_QC                              } from '../../modules/nf-core/biscuit/qc/main'
include { BISCUIT_VCF2BED as BISCUIT_VCF2BED_METH } from '../../modules/nf-core/biscuit/vcf2bed/main'
include { BISCUIT_MERGECG as BISCUIT_MERGECG_METH } from '../../modules/nf-core/biscuit/mergecg/main'
if (params.nomeseq) {
    include { BISCUIT_VCF2BED as BISCUIT_VCF2BED_NOME } from '../../modules/nf-core/biscuit/vcf2bed/main'
    include { BISCUIT_MERGECG as BISCUIT_MERGECG_NOME } from '../../modules/nf-core/biscuit/mergecg/main'
}

workflow BISCUIT {
    take:
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    versions = Channel.empty()

    /*
     * Generate biscuit index if not supplied
     */
    if (params.biscuit_index) {
        biscuit_index = file(params.biscuit_index)
    } else {
        BISCUIT_INDEX(params.fasta)
        biscuit_index = BISCUIT_INDEX.out.index
        versions = versions.mix(BISCUIT_INDEX.out.versions)
    }

    /*
     *  If we're using the --pbat option, biscuit will be set to -b 1.
     *  However, that would force the reads to align to the opposite strands we want,
     *  so we have to reverse the read order in the list (read 2 first, then read 1)
     */
    if (params.pbat) {
        reads = reads.map{meta, reads -> meta.single_end ? [meta, reads] : [meta, reads.reverse()]}
    }

    /*
     * Align with biscuit; mark duplicates unless params.skip_deduplication
     */
    if (params.skip_deduplication || params.rrbs){
        BISCUIT_ALIGN (
            reads,
            biscuit_index
        )
        versions = versions.mix(BISCUIT_ALIGN.out.versions)
        BISCUIT_ALIGN.out.bam
            .mix(BISCUIT_ALIGN.out.bai)
            .groupTuple(by: 0, size: 2, sort: {a, b -> a.toString() =~ /\.bai$/ ? 1: -1})
            .map{ meta, bam_bai -> [ meta, bam_bai[0], bam_bai[1] ] }
            .set{ alignments }
    } else {
        BISCUIT_BLASTER (
            reads,
            biscuit_index
        )
        versions = versions.mix(BISCUIT_BLASTER.out.versions)
        BISCUIT_BLASTER.out.bam
            .mix(BISCUIT_BLASTER.out.bai)
            .groupTuple(by: 0, size: 2, sort: {a, b -> a.toString() =~ /\.bai$/ ? 1: -1})
            .map{ meta, bam_bai -> [ meta, bam_bai[0], bam_bai[1] ] }
            .set{ alignments }
    }

    /*
     * Filter out reads with poor bisulfite conversion
     */
    if (params.bs_conv_filter < 1) {
        BISCUIT_BSCONV (
            alignments,
            biscuit_index
        )

        SAMTOOLS_INDEX_BSCONV ( BISCUIT_BSCONV.out.bsconv_bam )

        BISCUIT_BSCONV.out.bsconv_bam
            .mix(SAMTOOLS_INDEX_BSCONV.out.bai)
            .groupTuple(by: 0, size: 2, sort: {a, b -> a.toString() =~ /\.bai$/ ? 1: -1})
            .map{ meta, bam_bai -> [ meta, bam_bai[0], bam_bai[1] ] }
            .set{ alignments }
        versions = versions.mix(BISCUIT_BSCONV.out.versions)
        versions = versions.mix(SAMTOOLS_INDEX_BSCONV.out.versions)
    }

    /*
     * Extract all snp and methlyation information
     */
    BISCUIT_PILEUP (
        alignments.map{ meta, bam, bai -> [ meta, bam, bai, [], [] ] }, // add in blank lists for paired 'tumor' bam and index
        biscuit_index
    )
    versions = versions.mix(BISCUIT_PILEUP.out.versions)

    /*
     * Extract methlation information
     */
    BISCUIT_VCF2BED_METH (
        BISCUIT_PILEUP.out.vcf
    )
    versions = versions.mix(BISCUIT_VCF2BED_METH.out.versions)

    /*
    * Extract accessibility information
    */
    if (params.nomeseq) {
        BISCUIT_VCF2BED_NOME (
            BISCUIT_PILEUP.out.vcf
        )
        nome_bed = BISCUIT_VCF2BED_NOME.out.bed
        versions = versions.mix(BISCUIT_VCF2BED_NOME.out.versions)
    }

    /*
     * Merge information from neighboring C and G in CpG context
     */
    if (params.merge_cg) {
        BISCUIT_MERGECG_METH (
            BISCUIT_VCF2BED_METH.out.bed,
            biscuit_index
        )
        me_bed = BISCUIT_MERGECG_METH.out.mergecg_bed
        versions = versions.mix(BISCUIT_MERGECG_METH.out.versions)

        if (params.nomeseq) {
            BISCUIT_MERGECG_NOME (
                BISCUIT_VCF2BED_NOME.out.bed,
                biscuit_index
            )
            nome_bed = BISCUIT_MERGECG_NOME.out.mergecg_bed
            versions = versions.mix(BISCUIT_MERGECG_NOME.out.versions)
        }
    } else {
        me_bed = BISCUIT_VCF2BED_METH.out.bed
    }

    /*
     * QC alignments
     */
    BISCUIT_QC (
        alignments.map{meta, bam, bai -> [meta, bam]},
        biscuit_index
    )
    versions = versions.mix(BISCUIT_QC.out.versions)

    /*
     * Collect MultiQC inputs
     */
    if (!params.skip_multiqc) {
            BISCUIT_QC.out.biscuit_qc_reports.collect{ it[1] }
            .set{ multiqc_files }
    } else {
        multiqc_files = Channel.empty()
    }

    emit:
    bam   = alignments.map{meta, bam, bai -> [ meta, bam]}
    dedup = alignments.map{meta, bam, bai -> [ meta, bam]}
    mqc   = multiqc_files                        // path: *{html,txt}
    versions
}


// TODO:
// methylation extraction
// samtools flagstat, stats
// reports
// multiqc




