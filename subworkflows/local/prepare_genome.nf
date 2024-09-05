//
// Prepare reference genome files
//

include { UNTAR                       } from '../../modules/nf-core/untar/main'
include { GUNZIP                      } from '../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION   } from '../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX               } from '../../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX              } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    fasta         //      file: /path/to/genome.fasta
    fasta_index   //      file: /path/to/genome.fasta.fai
    bismark_index // directory: /path/to/bismark/index/
    bwameth_index // directory: /path/to/bwameth/index/

    main:
    ch_versions      = Channel.empty()
    ch_fasta         = Channel.empty()
    ch_fasta_index   = Channel.empty()
    ch_bismark_index = Channel.empty()
    ch_bwameth_index = Channel.empty()

    // FASTA, if supplied
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }
    ch_fasta.dump(tag: 'PREPARE_GENOME: ch_fasta')

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){

        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            if (bismark_index.endsWith('.gz')) {
                ch_bismark_index = UNTAR ( [ [:], file(bismark_index) ] ).untar.map { it[1] }
            } else {
                ch_bismark_index = Channel.value(file(bismark_index))
            }
        } else {
            BISMARK_GENOMEPREPARATION(ch_fasta)
            ch_bismark_index = BISMARK_GENOMEPREPARATION.out.index
            ch_versions = ch_versions.mix(BISMARK_GENOMEPREPARATION.out.versions)
        }
        ch_bismark_index.dump(tag: 'PREPARE_GENOME: ch_bismark_index')

    }
    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){

        /*
         * Generate bwameth index if not supplied
         */
        if (bwameth_index) {
            if (bwameth_index.endsWith('.tar.gz')) {
                ch_bismark_index = UNTAR ( [ [:], file(bwameth_index) ] ).untar.map { it[1] }
            } else {
                ch_bismark_index = Channel.value(file(bwameth_index))
            }
        } else {
            BWAMETH_INDEX(ch_fasta)
            ch_bwameth_index = BWAMETH_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMETH_INDEX.out.versions)
        }
        ch_bwameth_index.dump(tag: 'PREPARE_GENOME: ch_bwameth_index')

        /*
         * Generate fasta index if not supplied
         */
        if (fasta_index) {
            ch_fasta_index = Channel.value(file(fasta_index))
        } else {
            SAMTOOLS_FAIDX(ch_fasta.map{ fasta -> [[:], fasta]})
            ch_fasta_index = SAMTOOLS_FAIDX.out.fai.map{ return(it[1])}
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
        ch_fasta_index.dump(tag: 'PREPARE_GENOME: ch_fasta_index')
    }

    emit:
    fasta         = ch_fasta                  // channel: path(genome.fasta)
    bismark_index = ch_bismark_index          // channel: path(genome.fasta)
    bwameth_index = ch_bwameth_index          // channel: path(genome.fasta)
    fasta_index   = ch_fasta_index            // channel: path(genome.fasta)
    versions      = ch_versions               // channel: [ versions.yml ]

}
