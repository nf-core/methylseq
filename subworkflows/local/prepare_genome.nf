//
// Prepare reference genome files
//

include { BISMARK_GENOMEPREPARATION   } from '../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX               } from '../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list   : tools to prepare indices for
    biotype              // string : if additional fasta file is provided biotype value to use when appending entries to GTF file
    is_aws_igenome       // boolean: whether the genome files are from AWS iGenomes

    main:
    versions = Channel.empty()

    fasta = Channel.empty()
    bismark_index = Channel.empty()
    bwameth_index = Channel.empty()
    fasta_index = Channel.empty()

    // FASTA, if supplied
    if (params.fasta) {
        fasta = file(params.fasta)
    }

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){

        /*
         * Generate bismark index if not supplied
         */
        if (params.bismark_index) {
            bismark_index = file(params.bismark_index)
        } else {
            BISMARK_GENOMEPREPARATION(fasta)
            bismark_index = BISMARK_GENOMEPREPARATION.out.index
            versions = versions.mix(BISMARK_GENOMEPREPARATION.out.versions)
        }

    }
    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){

        /*
         * Generate bwameth index if not supplied
         */
        if (params.bwa_meth_index) {
            bwameth_index = file(params.bwa_meth_index)
        } else {
            BWAMETH_INDEX(fasta)
            bwameth_index = BWAMETH_INDEX.out.index
            versions = versions.mix(BWAMETH_INDEX.out.versions)
        }

        /*
         * Generate fasta index if not supplied
         */
        if (params.fasta_index) {
            fasta_index = file(params.fasta_index)
        } else {
            SAMTOOLS_FAIDX([[:], fasta])
            fasta_index = SAMTOOLS_FAIDX.out.fai.map{ return(it[1])}
            versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    }

    emit:
    fasta = fasta
    bismark_index = bismark_index
    bwameth_index = bwameth_index
    fasta_index = fasta_index
    versions = versions

}
