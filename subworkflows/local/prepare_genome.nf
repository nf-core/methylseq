//
// Prepare reference genome files
//

include { UNTAR                       } from '../../modules/nf-core/untar/main'
include { BISMARK_GENOMEPREPARATION   } from '../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX               } from '../../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX              } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    main:
    ch_versions      = Channel.empty()
    ch_fasta         = Channel.empty()
    ch_bismark_index = Channel.empty()
    ch_bwameth_index = Channel.empty()
    ch_fasta_index   = Channel.empty()

    // FASTA, if supplied
    if (params.fasta) {
        ch_fasta = Channel.value(file(params.fasta))
    }

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){

        /*
         * Generate bismark index if not supplied
         */
        if (params.bismark_index) {
            if (params.bismark_index.endsWith('.gz')) {
                ch_bismark_index = UNTAR ( [ [:], file(params.bismark_index) ] ).untar.map { it[1] }
            } else {
                ch_bismark_index = Channel.value(file(params.bismark_index))
            }
        } else {
            BISMARK_GENOMEPREPARATION(ch_fasta)
            ch_bismark_index = BISMARK_GENOMEPREPARATION.out.index
            ch_versions = ch_versions.mix(BISMARK_GENOMEPREPARATION.out.versions)
        }

    }
    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){

        /*
         * Generate bwameth index if not supplied
         */
        if (params.bwa_meth_index) {
            if (params.bwa_meth_index.endsWith('.tar.gz')) {
                ch_bismark_index = UNTAR ( [ [:], file(params.bwa_meth_index) ] ).untar.map { it[1] }
            } else {
                ch_bismark_index = Channel.value(file(params.bwa_meth_index))
            }
        } else {
            BWAMETH_INDEX(ch_fasta)
            ch_bwameth_index = BWAMETH_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMETH_INDEX.out.versions)
        }

        /*
         * Generate fasta index if not supplied
         */
        if (params.fasta_index) {
            ch_fasta_index = Channel.value(file(params.fasta_index))
        } else {
            SAMTOOLS_FAIDX([[:], ch_fasta])
            ch_fasta_index = SAMTOOLS_FAIDX.out.fai.map{ return(it[1])}
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta                  // channel: path(genome.fasta)
    bismark_index = ch_bismark_index          // channel: path(genome.fasta)
    bwameth_index = ch_bwameth_index          // channel: path(genome.fasta)
    fasta_index   = ch_fasta_index            // channel: path(genome.fasta)
    versions      = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}
