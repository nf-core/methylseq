/*
 * biscuit subworkflow
 */

include { BISCUIT_INDEX   } from '../../modules/nf-core/biscuit/index/main'
include { BISCUIT_ALIGN   } from '../../modules/nf-core/biscuit/align/main'
include { BISCUIT_BLASTER } from '../../modules/nf-core/biscuit/biscuitblaster/main'
include { BISCUIT_BSCONV  } from '../../modules/nf-core/biscuit/bsconv/main'
include { BISCUIT_PILEUP  } from '../../modules/nf-core/biscuit/pileup/main'
include { BISCUIT_VCF2BED } from '../../modules/nf-core/biscuit/vcf2bed/main'
include { BISCUIT_MERGECG } from '../../modules/nf-core/biscuit/mergecg/main'
include { BISCUIT_QC      } from '../../modules/nf-core/biscuit/qc/main'
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
        alignments = BISCUIT_ALIGN.out.bam
    } else {
        BISCUIT_BLASTER (
            reads,
            biscuit_index
        )
        versions = versions.mix(BISCUIT_BLASTER.out.versions)
        alignments = BISCUIT_BLASTER.out.bam
    }

    /*
     * Filter out reads with poor bisulfite conversion
     */
    if (params.skip_bs_conv_filter < 1) {
        BISCUIT_BSCONV(
            alignments,
            biscuit_index
        )
        alignments = BISCUIT_BSCONV.out.bsconv_bam
        versions = versions.mix(BISCUIT_BSCONV.out.versions)
    }

    /*
     * Extract all snp and methlyation information
     */
    BISCUIT_PILEUP (
        alignments,
        biscuit_index
    )
    versions = versions.mix(BISCUIT_PILEUP.out.versions)

    /*
     * Extract methlation information
     */
    BISCUIT_VCF2BED (
        BISCUIT_PILEUP.out.vcf
    )
    versions = versions.mix(BISCUIT_VCF2BED.out.versions)

    if (params.nomeseq) {
        /*
        * Extract accessibility information
        */
        BISCUIT_VCF2BED_NOME (
            BISCUIT_PILEUP.out.vcf
        )
        versions = versions.mix(BISCUIT_VCF2BED_NOME.out.versions)

        /*
        * Merge information from neighboring C and G, where both are in HCGD context
        */
        if (params.merge_cg) {
            BISCUIT_MERGECG_NOME (
                BISCUIT_VCF2BED.out.bed,
                biscuit_index
            )
            me_bed = BISCUIT_MERGECG.out.bed
            versions = versions.mix(BISCUIT_MERGECG.out.versions)
        } else {
            me_bed = BISCUIT_VCF2BED.out.bed
        }

    }

    /*
     * Merge information from neighboring C and G in CpG context
     */
    if (params.merge_cg) {
        BISCUIT_MERGECG (
            BISCUIT_VCF2BED.out.bed,
            biscuit_index
        )
        me_bed = BISCUIT_MERGECG.out.bed
        versions = versions.mix(BISCUIT_MERGECG.out.versions)
    } else {
        me_bed = BISCUIT_VCF2BED.out.bed
    }

    /*
     * QC alignments
     */
    BISCUIT_QC (
        alignments,
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




